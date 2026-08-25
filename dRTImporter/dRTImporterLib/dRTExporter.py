"""Fast dynamic RT Structure Set exporter for SlicerDynamicPET.

This exporter runs inside 3D Slicer. It reads the binary labelmap representation
of each temporal segment directly, extracts planar contours with VTK marching
squares, converts Slicer RAS coordinates to DICOM LPS coordinates, and writes
one temporal RTSTRUCT with pydicom.

Unlike the earlier implementation, export does NOT round-trip a large temporary
segmentation through SlicerRT. SlicerRT remains required by dRTImporter for
RTSTRUCT import, but is not on the export hot path.
"""

import faulthandler
import hashlib
import json
import logging
import math
import os
import sys
import struct
import tempfile
import time
from datetime import datetime


ENHANCED_PET_IMAGE_STORAGE_UID = "1.2.840.10008.5.1.4.1.1.130"
RTSTRUCT_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.481.3"
RT_REFERENCED_STUDY_SOP_CLASS_UID = "1.2.840.10008.3.1.2.3.2"


_DEBUG_LOG_PATH = os.path.join(tempfile.gettempdir(), 'dRTExporter_debug.log')
_FAULT_LOG_PATH = os.path.join(tempfile.gettempdir(), 'dRTExporter_crash.log')
_debug_stream = None
_fault_stream = None
_debug_start_time = None
_fault_handler_was_enabled = False
_DIAGNOSTICS_ENABLED = str(os.environ.get('SLICER_DRT_DEBUG', '0')).lower() in ('1', 'true', 'yes', 'on')


def _start_diagnostics():
  """Start line-buffered crash diagnostics when SLICER_DRT_DEBUG=1."""
  if not _DIAGNOSTICS_ENABLED:
    return
  global _debug_stream, _fault_stream, _debug_start_time, _fault_handler_was_enabled
  _debug_start_time = time.perf_counter()
  _fault_handler_was_enabled = faulthandler.is_enabled()
  try:
    _debug_stream = open(_DEBUG_LOG_PATH, 'w', buffering=1)
  except Exception:
    _debug_stream = None
  try:
    _fault_stream = open(_FAULT_LOG_PATH, 'w', buffering=1)
    faulthandler.enable(file=_fault_stream, all_threads=True)
  except Exception:
    _fault_stream = None
  _dbg('diagnostics started')
  _dbg(f'debug log: {_DEBUG_LOG_PATH}')
  _dbg(f'fault log: {_FAULT_LOG_PATH}')


def _stop_diagnostics():
  if not _DIAGNOSTICS_ENABLED:
    return
  global _debug_stream, _fault_stream, _fault_handler_was_enabled
  _dbg('diagnostics stopped')
  try:
    if _debug_stream is not None:
      _debug_stream.flush()
      _debug_stream.close()
  except Exception:
    pass
  _debug_stream = None
  try:
    if _fault_stream is not None:
      faulthandler.disable()
      _fault_stream.flush()
      _fault_stream.close()
      if _fault_handler_was_enabled:
        faulthandler.enable(file=sys.stderr, all_threads=True)
  except Exception:
    pass
  _fault_stream = None


def _dbg(message):
  if not _DIAGNOSTICS_ENABLED:
    return
  elapsed = (time.perf_counter() - _debug_start_time) if _debug_start_time else 0.0
  text = f'[dRT DEBUG +{elapsed:8.3f}s] {message}'
  try:
    print(text, flush=True)
  except Exception:
    pass
  try:
    if _debug_stream is not None:
      _debug_stream.write(text + '\n')
      _debug_stream.flush()
  except Exception:
    pass


def _segment_ids(segmentation_node):
  import vtk

  ids = vtk.vtkStringArray()
  segmentation_node.GetSegmentation().GetSegmentIDs(ids)
  return [ids.GetValue(index) for index in range(ids.GetNumberOfValues())]


def _parse_dicom_vector(value, expected_length=None):
  if value in (None, ''):
    return None
  try:
    if isinstance(value, (list, tuple)):
      values = [float(item) for item in value]
    else:
      values = [float(item) for item in str(value).split('\\') if item != '']
    if expected_length is not None and len(values) != expected_length:
      return None
    return values
  except Exception:
    return None


def _validate_pet_reference_metadata(metadata, expected_count, require_enhanced_pet=False):
  if not isinstance(metadata, dict):
    raise ValueError('PET DICOM provenance is not a dictionary.')
  frames = metadata.get('frames') or []
  if len(frames) != expected_count:
    raise ValueError(
      f'PET DICOM provenance has {len(frames)} frames; expected {expected_count}.')
  for keyword in ('studyInstanceUID', 'seriesInstanceUID', 'frameOfReferenceUID'):
    if not str(metadata.get(keyword) or ''):
      raise ValueError(f'PET DICOM provenance is missing {keyword}.')

  seen_uids = set()
  for frame_index, frame in enumerate(frames):
    instances = frame.get('instances') or []
    if not instances:
      raise ValueError(f'PET temporal frame {frame_index + 1} has no source SOP instances.')
    for instance in instances:
      uid = str(instance.get('sopInstanceUID') or '')
      sop_class = str(instance.get('sopClassUID') or '')
      if not uid or not sop_class:
        raise ValueError(
          f'PET temporal frame {frame_index + 1} has incomplete SOP identity.')
      if uid in seen_uids:
        raise ValueError(
          'The same PET SOP Instance UID occurs in more than one temporal frame.')
      seen_uids.add(uid)

  if require_enhanced_pet:
    for frame_index, frame in enumerate(frames):
      instances = frame.get('instances') or []
      if (len(instances) != 1
          or str(instances[0].get('sopClassUID') or '')
          != ENHANCED_PET_IMAGE_STORAGE_UID):
        raise ValueError(
          'require_enhanced_pet=True requires exactly one Enhanced PET SOP '
          f'instance per temporal frame; frame {frame_index + 1} does not match.')
  return metadata


def _metadata_from_dpet_importer(pet_sequence_node):
  text = pet_sequence_node.GetAttribute('dPET.DICOM.FrameReferences')
  if not text:
    return None
  try:
    metadata = json.loads(text)
    _validate_pet_reference_metadata(
      metadata, pet_sequence_node.GetNumberOfDataNodes(), False)
    _dbg('PET provenance source=dPETImporter persisted metadata')
    return metadata
  except Exception as error:
    _dbg(f'dPETImporter provenance present but invalid: {error}')
    return None


def _metadata_from_mrml_and_dicom_database(pet_sequence_node):
  """Recover provenance without reopening DICOM files when MRML/DB is sufficient."""
  import slicer

  db = getattr(slicer, 'dicomDatabase', None)
  frame_count = pet_sequence_node.GetNumberOfDataNodes()
  frames = []

  def db_value(uid, tag):
    if db is None or not uid:
      return ''
    try:
      file_path = db.fileForInstance(uid)
      if not file_path:
        return ''
      value = db.fileValue(file_path, tag)
      return str(value) if value not in (None, '') else ''
    except Exception:
      return ''

  all_uids = []
  for frame_index in range(frame_count):
    node = pet_sequence_node.GetNthDataNode(frame_index)
    uids = (node.GetAttribute('DICOM.instanceUIDs') or '').split() if node else []
    if not uids:
      return None
    instances = []
    for uid in uids:
      sop_class = (node.GetAttribute('DICOM.SOPClassUID') or '') if len(uids) == 1 else ''
      sop_class = sop_class or db_value(uid, '0008,0016')
      if not sop_class:
        return None
      item = {
        'sopInstanceUID': uid,
        'sopClassUID': sop_class,
      }
      position = _parse_dicom_vector(db_value(uid, '0020,0032'), 3)
      orientation = _parse_dicom_vector(db_value(uid, '0020,0037'), 6)
      if position is not None:
        item['imagePositionPatient'] = position
      if orientation is not None:
        item['imageOrientationPatient'] = orientation
      instances.append(item)
      all_uids.append(uid)
    frames.append({'index': frame_index, 'instances': instances})

  first_node = pet_sequence_node.GetNthDataNode(0) if frame_count else None
  first_uid = all_uids[0] if all_uids else ''
  study_uid = (
    pet_sequence_node.GetAttribute('dPET.DICOM.StudyInstanceUID')
    or pet_sequence_node.GetAttribute('DICOM.StudyInstanceUID')
    or (first_node.GetAttribute('DICOM.StudyInstanceUID') if first_node else '')
    or db_value(first_uid, '0020,000D'))
  series_uid = (
    pet_sequence_node.GetAttribute('dPET.DICOM.SeriesInstanceUID')
    or pet_sequence_node.GetAttribute('dPET.SourceSeriesInstanceUID')
    or pet_sequence_node.GetAttribute('DICOM.SeriesInstanceUID')
    or (first_node.GetAttribute('DICOM.SeriesInstanceUID') if first_node else '')
    or db_value(first_uid, '0020,000E'))
  frame_uid = (
    pet_sequence_node.GetAttribute('dPET.DICOM.FrameOfReferenceUID')
    or pet_sequence_node.GetAttribute('DICOM.FrameOfReferenceUID')
    or (first_node.GetAttribute('DICOM.FrameOfReferenceUID') if first_node else '')
    or db_value(first_uid, '0020,0052'))
  if not study_uid or not series_uid or not frame_uid:
    return None

  patient_study = {}
  for keyword, tag in {
      'SpecificCharacterSet': '0008,0005',
      'PatientName': '0010,0010', 'PatientID': '0010,0020',
      'PatientBirthDate': '0010,0030', 'PatientSex': '0010,0040',
      'PatientAge': '0010,1010', 'PatientSize': '0010,1020',
      'PatientWeight': '0010,1030', 'StudyDate': '0008,0020',
      'StudyTime': '0008,0030', 'AccessionNumber': '0008,0050',
      'StudyID': '0020,0010', 'ReferringPhysicianName': '0008,0090',
      'InstitutionName': '0008,0080', 'InstitutionAddress': '0008,0081',
      }.items():
    value = db_value(first_uid, tag)
    if value:
      patient_study[keyword] = value

  metadata = {
    'schemaVersion': 1,
    'studyInstanceUID': str(study_uid),
    'seriesInstanceUID': str(series_uid),
    'frameOfReferenceUID': str(frame_uid),
    'patientStudy': patient_study,
    'frames': frames,
  }
  try:
    _validate_pet_reference_metadata(metadata, frame_count, False)
  except Exception:
    return None
  _dbg('PET provenance source=MRML + DICOM database metadata')
  return metadata


def _metadata_from_source_files(pet_sequence_node, explicit_files=None):
  """Last-resort header-only fallback for sequences not loaded by dPETImporter."""
  import pydicom

  frame_count = pet_sequence_node.GetNumberOfDataNodes()
  frame_files = None
  if explicit_files is not None:
    supplied = list(explicit_files)
    if supplied and isinstance(supplied[0], (list, tuple)):
      frame_files = [[os.path.abspath(str(path)) for path in group] for group in supplied]
    elif len(supplied) == frame_count:
      frame_files = [[os.path.abspath(str(path))] for path in supplied]
    else:
      raise ValueError(
        'pet_dicom_files must be one file per temporal frame or a nested '
        'sequence of source files for each temporal frame.')
  else:
    flat_text = pet_sequence_node.GetAttribute('MultiVolume.FrameFileList')
    if not flat_text:
      return None
    flat = [path for path in flat_text.split(',') if path]
    if not flat or len(flat) % frame_count != 0:
      return None
    files_per_frame = len(flat) // frame_count
    frame_files = [
      flat[index * files_per_frame:(index + 1) * files_per_frame]
      for index in range(frame_count)
    ]

  first_dataset = None
  frames = []
  for frame_index, paths in enumerate(frame_files):
    instances = []
    for path in paths:
      if not os.path.isfile(path):
        return None
      dataset = pydicom.dcmread(path, stop_before_pixels=True, force=True)
      if first_dataset is None:
        first_dataset = dataset
      item = {
        'sopInstanceUID': str(getattr(dataset, 'SOPInstanceUID', '') or ''),
        'sopClassUID': str(getattr(dataset, 'SOPClassUID', '') or ''),
      }
      if hasattr(dataset, 'ImagePositionPatient'):
        item['imagePositionPatient'] = [float(value) for value in dataset.ImagePositionPatient]
      if hasattr(dataset, 'ImageOrientationPatient'):
        item['imageOrientationPatient'] = [float(value) for value in dataset.ImageOrientationPatient]
      instances.append(item)
    frames.append({'index': frame_index, 'instances': instances})

  if first_dataset is None:
    return None
  patient_study = {}
  for keyword in (
      'SpecificCharacterSet', 'PatientName', 'PatientID', 'PatientBirthDate',
      'PatientSex', 'PatientAge', 'PatientSize', 'PatientWeight', 'StudyDate',
      'StudyTime', 'AccessionNumber', 'StudyID', 'ReferringPhysicianName',
      'InstitutionName', 'InstitutionAddress'):
    if hasattr(first_dataset, keyword):
      patient_study[keyword] = str(getattr(first_dataset, keyword))

  metadata = {
    'schemaVersion': 1,
    'studyInstanceUID': str(getattr(first_dataset, 'StudyInstanceUID', '') or ''),
    'seriesInstanceUID': str(getattr(first_dataset, 'SeriesInstanceUID', '') or ''),
    'frameOfReferenceUID': str(getattr(first_dataset, 'FrameOfReferenceUID', '') or ''),
    'patientStudy': patient_study,
    'frames': frames,
  }
  _validate_pet_reference_metadata(metadata, frame_count, False)
  _dbg('PET provenance source=source DICOM header fallback')
  return metadata


def _resolve_pet_reference_metadata(
    pet_sequence_node, explicit_files=None, require_enhanced_pet=False):
  frame_count = pet_sequence_node.GetNumberOfDataNodes()
  metadata = _metadata_from_dpet_importer(pet_sequence_node)
  if metadata is None:
    metadata = _metadata_from_mrml_and_dicom_database(pet_sequence_node)
  if metadata is None:
    metadata = _metadata_from_source_files(pet_sequence_node, explicit_files)
  if metadata is None:
    raise ValueError(
      'Unable to resolve PET DICOM provenance. Export works best with '
      'dPETImporter metadata, but can also use MRML/DICOM database metadata '
      'or source DICOM headers when available.')

  _validate_pet_reference_metadata(metadata, frame_count, require_enhanced_pet)

  # Cache successful fallback resolution in the scene so subsequent exports and
  # MRML scene saves no longer depend on the original source files.
  if not pet_sequence_node.GetAttribute('dPET.DICOM.FrameReferences'):
    pet_sequence_node.SetAttribute(
      'dPET.DICOM.FrameReferences',
      json.dumps(metadata, separators=(',', ':'), ensure_ascii=False))
    pet_sequence_node.SetAttribute(
      'dPET.DICOM.StudyInstanceUID', metadata['studyInstanceUID'])
    pet_sequence_node.SetAttribute(
      'dPET.DICOM.SeriesInstanceUID', metadata['seriesInstanceUID'])
    pet_sequence_node.SetAttribute(
      'dPET.DICOM.FrameOfReferenceUID', metadata['frameOfReferenceUID'])
  return metadata


def _reference_frame_of_reference_uid(reference_volume_node):
  """Resolve reference-volume Frame of Reference without reopening DICOM."""
  import slicer

  if reference_volume_node is None:
    return None
  direct = reference_volume_node.GetAttribute('DICOM.FrameOfReferenceUID')
  if direct:
    return str(direct)
  instance_uids = (reference_volume_node.GetAttribute('DICOM.instanceUIDs') or '').split()
  if not instance_uids:
    return None
  try:
    file_path = slicer.dicomDatabase.fileForInstance(instance_uids[0])
    if not file_path:
      return None
    value = slicer.dicomDatabase.fileValue(file_path, '0020,0052')
    return str(value) if value not in (None, '') else None
  except Exception:
    return None


def _segment_definitions(segmentation_sequence_node):
  """Return stable segment definitions keyed by persistent Slicer segment ID."""
  definitions = {}
  for item_index in range(segmentation_sequence_node.GetNumberOfDataNodes()):
    segmentation_node = segmentation_sequence_node.GetNthDataNode(item_index)
    for segment_id in _segment_ids(segmentation_node):
      segment = segmentation_node.GetSegmentation().GetSegment(segment_id)
      if segment is None:
        continue
      if segment_id not in definitions:
        definitions[segment_id] = {
          'id': segment_id,
          'name': segment.GetName(),
          'color': tuple(float(value) for value in segment.GetColor()),
        }
  return definitions


def _binary_labelmap(segmentation_node, segment_id):
  """Get one isolated binary labelmap in world geometry."""
  import slicer
  import vtkSegmentationCorePython as vtkSegmentationCore

  _dbg(f'_binary_labelmap BEGIN segment={segment_id}')
  labelmap = vtkSegmentationCore.vtkOrientedImageData()
  binary_name = (
    vtkSegmentationCore.vtkSegmentationConverter
    .GetSegmentationBinaryLabelmapRepresentationName())
  segment = segmentation_node.GetSegmentation().GetSegment(segment_id)

  # Prefer the MRML node API when the binary representation already exists.
  # It explicitly returns a copy containing only the requested segment, even
  # when several segments share one internal binary-labelmap layer.
  if segment is not None and segment.GetRepresentation(binary_name) is not None:
    _dbg('vtkMRMLSegmentationNode.GetBinaryLabelmapRepresentation BEGIN')
    success = segmentation_node.GetBinaryLabelmapRepresentation(
      segment_id, labelmap)
    _dbg(
      'vtkMRMLSegmentationNode.GetBinaryLabelmapRepresentation END '
      f'success={bool(success)}')
  else:
    # If binary labelmap is not present (for example a Planar-contour source),
    # use Slicer logic to create/extract only this requested segment.
    _dbg('vtkSlicerSegmentationsModuleLogic.GetSegmentBinaryLabelmapRepresentation BEGIN')
    success = slicer.vtkSlicerSegmentationsModuleLogic.GetSegmentBinaryLabelmapRepresentation(
      segmentation_node, segment_id, labelmap)
    _dbg(
      'vtkSlicerSegmentationsModuleLogic.GetSegmentBinaryLabelmapRepresentation END '
      f'success={bool(success)}')
  if not success:
    return None
  _dbg(f'labelmap extent={tuple(labelmap.GetExtent())} dims={tuple(labelmap.GetDimensions())}')
  scalars = labelmap.GetPointData().GetScalars()
  if scalars is None or scalars.GetNumberOfTuples() == 0:
    _dbg('_binary_labelmap END empty scalars')
    return None
  scalar_range = scalars.GetRange()
  _dbg(f'labelmap scalarType={scalars.GetDataTypeAsString()} range={scalar_range}')
  if scalar_range[1] <= 0.0:
    _dbg('_binary_labelmap END zero mask')
    return None
  _dbg('_binary_labelmap END non-empty')
  return labelmap


def _matrix_values(matrix):
  return tuple(
    round(float(matrix.GetElement(row, column)), 12)
    for row in range(4)
    for column in range(4)
  )


def _labelmap_signature(labelmap):
  """Exact content+geometry fingerprint used to reuse unchanged contours."""
  import vtk
  from vtk.util.numpy_support import vtk_to_numpy

  _dbg('_labelmap_signature BEGIN')
  scalars = labelmap.GetPointData().GetScalars()
  if scalars is None:
    _dbg('_labelmap_signature END no scalars')
    return None

  matrix = vtk.vtkMatrix4x4()
  labelmap.GetImageToWorldMatrix(matrix)
  extent = tuple(int(value) for value in labelmap.GetExtent())
  _dbg(f'_labelmap_signature vtk_to_numpy BEGIN extent={extent}')
  array = vtk_to_numpy(scalars)
  _dbg(f'_labelmap_signature vtk_to_numpy END size={array.size} bytes={array.nbytes}')

  digest = hashlib.blake2b(digest_size=16)
  digest.update(str(extent).encode('ascii'))
  digest.update(repr(_matrix_values(matrix)).encode('ascii'))
  digest.update(memoryview(array).cast('B'))
  result = digest.hexdigest()
  _dbg(f'_labelmap_signature END signature={result}')
  return result


def _remove_duplicate_and_collinear_points(points, tolerance=1e-7):
  """Reduce a closed contour without changing its polygonal geometry."""
  import math

  if len(points) < 3:
    return points

  # Remove an explicit repeated closing point. CLOSED_PLANAR implicitly closes
  # the polygon from the last point back to the first.
  def distance2(a, b):
    return sum((a[i] - b[i]) ** 2 for i in range(3))

  cleaned = []
  for point in points:
    if not cleaned or distance2(point, cleaned[-1]) > tolerance * tolerance:
      cleaned.append(point)
  if len(cleaned) >= 2 and distance2(cleaned[0], cleaned[-1]) <= tolerance * tolerance:
    cleaned.pop()

  if len(cleaned) < 4:
    return cleaned

  # Marching squares produces long exactly-collinear runs along voxel edges.
  # Removing middle vertices is lossless and keeps ContourData much smaller.
  changed = True
  while changed and len(cleaned) >= 4:
    changed = False
    keep = []
    count = len(cleaned)
    for index in range(count):
      previous = cleaned[(index - 1) % count]
      current = cleaned[index]
      following = cleaned[(index + 1) % count]
      a = [current[i] - previous[i] for i in range(3)]
      b = [following[i] - current[i] for i in range(3)]
      cross = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
      ]
      cross_norm = math.sqrt(sum(value * value for value in cross))
      scale = max(
        1.0,
        math.sqrt(sum(value * value for value in a)),
        math.sqrt(sum(value * value for value in b)),
      )
      if cross_norm <= tolerance * scale:
        changed = True
      else:
        keep.append(current)
    if len(keep) < 3:
      break
    cleaned = keep
  return cleaned


def _contours_from_binary_labelmap(labelmap):
  """Extract CLOSED_PLANAR contours with vtkMarchingSquares on explicit 2D slices.

  The previous diagnostic build drove one vtkMarchingSquares instance with a
  3D vtkImageData input and changed ImageRange for every slice.  A reproducible
  native crash was localized to marching.Update() in that 3D ImageRange path.

  This implementation deliberately keeps vtkMarchingSquares + vtkStripper as
  the contour algorithm, but gives MarchingSquares a true 2D vtkImageData for
  each non-empty slice.  Therefore SetImageRange is not used at all.
  """
  import numpy as np
  import vtk
  from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

  _dbg('_contours_from_binary_labelmap BEGIN explicit-2D mode')
  image_to_ras = vtk.vtkMatrix4x4()
  labelmap.GetImageToWorldMatrix(image_to_ras)

  scalars = labelmap.GetPointData().GetScalars()
  if scalars is None or scalars.GetNumberOfTuples() == 0:
    _dbg('_contours END no scalars')
    return []

  scalar_range = scalars.GetRange()
  if scalar_range[1] <= 0.0:
    _dbg('_contours END zero range')
    return []

  # Keep the same binary thresholding semantics as the diagnostic build.
  threshold = vtk.vtkImageThreshold()
  threshold.SetInputData(labelmap)
  threshold.ThresholdBetween(1.0, scalar_range[1])
  threshold.SetInValue(1)
  threshold.SetOutValue(0)
  threshold.ReplaceInOn()
  threshold.ReplaceOutOn()
  threshold.SetOutputScalarTypeToUnsignedChar()
  _dbg('vtkImageThreshold.Update BEGIN')
  threshold.Update()
  _dbg('vtkImageThreshold.Update END')

  binary_image = threshold.GetOutput()
  extent = binary_image.GetExtent()
  dimensions = binary_image.GetDimensions()
  _dbg(f'binary image extent={tuple(extent)} dims={tuple(dimensions)}')
  if dimensions[0] <= 0 or dimensions[1] <= 0 or dimensions[2] <= 0:
    _dbg('_contours END invalid dimensions')
    return []

  _dbg('binary vtk_to_numpy BEGIN')
  array = vtk_to_numpy(binary_image.GetPointData().GetScalars())
  _dbg(f'binary vtk_to_numpy END size={array.size}')
  array = array.reshape((dimensions[2], dimensions[1], dimensions[0]))
  nonzero_count = int(np.count_nonzero(array))
  nonempty_local_slices = np.flatnonzero(np.any(array != 0, axis=(1, 2)))
  _dbg(
    f'binary nonzero={nonzero_count} '
    f'nonemptyLocalSlices={nonempty_local_slices.tolist()}')
  if nonempty_local_slices.size == 0:
    _dbg('_contours END no nonempty slices')
    return []

  contours = []
  width = dimensions[0]
  height = dimensions[1]

  for local_k_value in nonempty_local_slices:
    local_k = int(local_k_value)
    original_k = extent[4] + local_k
    _dbg(
      f'slice localK={local_k} originalK={original_k} '
      f'explicit 2D preparation BEGIN')

    # Pad by one zero-valued pixel in X/Y so an ROI touching the cropped
    # labelmap boundary still produces a closed contour.  This replaces the
    # former 3D vtkImageConstantPad + SetImageRange combination.
    slice_array = np.ascontiguousarray(array[local_k], dtype=np.uint8)
    padded_slice = np.pad(
      slice_array,
      ((1, 1), (1, 1)),
      mode='constant',
      constant_values=0)

    image2d = vtk.vtkImageData()
    image2d.SetDimensions(width + 2, height + 2, 1)
    image2d.SetExtent(0, width + 1, 0, height + 1, 0, 0)
    image2d.SetOrigin(0.0, 0.0, 0.0)
    image2d.SetSpacing(1.0, 1.0, 1.0)
    vtk_scalars = numpy_to_vtk(
      num_array=padded_slice.ravel(order='C'),
      deep=True,
      array_type=vtk.VTK_UNSIGNED_CHAR)
    image2d.GetPointData().SetScalars(vtk_scalars)

    _dbg(
      f'slice originalK={original_k} explicit 2D preparation END '
      f'dims={image2d.GetDimensions()}')

    # Create fresh filters per slice.  This avoids carrying any internal
    # locator/pipeline state across successive temporal planes.
    marching = vtk.vtkMarchingSquares()
    marching.SetInputData(image2d)
    marching.SetValue(0, 0.5)
    _dbg(f'slice originalK={original_k} vtkMarchingSquares.Update BEGIN')
    marching.Update()
    marching_output = marching.GetOutput()
    _dbg(
      f'slice originalK={original_k} vtkMarchingSquares.Update END '
      f'points={marching_output.GetNumberOfPoints()} '
      f'cells={marching_output.GetNumberOfCells()}')

    if marching_output.GetNumberOfCells() == 0:
      continue

    stripper = vtk.vtkStripper()
    stripper.SetInputData(marching_output)
    stripper.JoinContiguousSegmentsOn()
    _dbg(f'slice originalK={original_k} vtkStripper.Update BEGIN')
    stripper.Update()
    polydata = stripper.GetOutput()
    _dbg(
      f'slice originalK={original_k} vtkStripper.Update END '
      f'points={polydata.GetNumberOfPoints()} '
      f'cells={polydata.GetNumberOfCells()}')

    _dbg(f'slice originalK={original_k} contour cell extraction BEGIN')
    for cell_index in range(polydata.GetNumberOfCells()):
      cell = polydata.GetCell(cell_index)
      if cell is None or cell.GetNumberOfPoints() < 3:
        continue

      points_lps = []
      for point_index in range(cell.GetNumberOfPoints()):
        point_id = cell.GetPointId(point_index)
        local_point = polydata.GetPoint(point_id)

        # MarchingSquares coordinates are in the padded local 2D image.
        # Convert them back to the original labelmap IJK before applying the
        # labelmap's ImageToWorld matrix.  The -1 removes the explicit pad.
        original_i = extent[0] + float(local_point[0]) - 1.0
        original_j = extent[2] + float(local_point[1]) - 1.0
        ras4 = image_to_ras.MultiplyPoint(
          (original_i, original_j, float(original_k), 1.0))
        points_lps.append(
          (-float(ras4[0]), -float(ras4[1]), float(ras4[2])))

      points_lps = _remove_duplicate_and_collinear_points(points_lps)
      if len(points_lps) >= 3:
        contours.append(points_lps)
    _dbg(
      f'slice originalK={original_k} contour cell extraction END '
      f'totalContours={len(contours)}')

  _dbg(f'_contours_from_binary_labelmap END contours={len(contours)}')
  return contours


def _copy_if_present(source, destination, keywords):
  for keyword in keywords:
    if hasattr(source, keyword):
      setattr(destination, keyword, getattr(source, keyword))


def _new_rtstruct_dataset(pet_metadata, series_description, series_number, structure_set_label):
  """Create a compact RTSTRUCT using persisted PET patient/study provenance."""
  from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
  from pydicom.sequence import Sequence
  from pydicom.uid import ExplicitVRLittleEndian, PYDICOM_IMPLEMENTATION_UID, generate_uid

  now = datetime.now()
  sop_instance_uid = generate_uid()
  file_meta = FileMetaDataset()
  file_meta.FileMetaInformationVersion = b'\x00\x01'
  file_meta.MediaStorageSOPClassUID = RTSTRUCT_SOP_CLASS_UID
  file_meta.MediaStorageSOPInstanceUID = sop_instance_uid
  file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
  file_meta.ImplementationClassUID = PYDICOM_IMPLEMENTATION_UID

  rtstruct = FileDataset('', {}, file_meta=file_meta, preamble=b'\0' * 128)
  patient_study = pet_metadata.get('patientStudy') or {}
  rtstruct.SpecificCharacterSet = patient_study.get('SpecificCharacterSet') or 'ISO_IR 192'
  rtstruct.SOPClassUID = RTSTRUCT_SOP_CLASS_UID
  rtstruct.SOPInstanceUID = sop_instance_uid
  rtstruct.Modality = 'RTSTRUCT'

  # RTSTRUCT IOD Type 2 attributes must be present even when unknown.  Keeping
  # them as empty strings avoids the repeated SlicerRT/DCMTK conformance
  # warnings seen with the previous exporter.
  type2_defaults = {
    'PatientName': '',
    'PatientID': '',
    'PatientBirthDate': '',
    'PatientSex': '',
    'StudyDate': '',
    'StudyTime': '',
    'ReferringPhysicianName': '',
    'StudyID': '',
    'AccessionNumber': '',
    'OperatorsName': '',
  }
  for keyword, default_value in type2_defaults.items():
    value = patient_study.get(keyword)
    setattr(rtstruct, keyword, default_value if value in (None, '') else value)

  # Optional patient/study fields are copied only when available.
  for keyword in (
      'PatientAge', 'PatientSize', 'PatientWeight', 'InstitutionName',
      'InstitutionAddress', 'PerformingPhysicianName'):
    value = patient_study.get(keyword)
    if value not in (None, ''):
      setattr(rtstruct, keyword, value)

  rtstruct.StudyInstanceUID = str(pet_metadata['studyInstanceUID'])
  rtstruct.SeriesInstanceUID = generate_uid()
  rtstruct.SeriesNumber = int(series_number)
  rtstruct.InstanceNumber = 1
  rtstruct.SeriesDescription = str(series_description)
  rtstruct.SeriesDate = now.strftime('%Y%m%d')
  rtstruct.SeriesTime = now.strftime('%H%M%S.%f')
  rtstruct.StructureSetLabel = str(structure_set_label)[:16]
  rtstruct.StructureSetName = str(series_description)
  rtstruct.StructureSetDescription = (
    'Dynamic RTSTRUCT generated by 3D Slicer / SlicerDynamicPET')
  rtstruct.StructureSetDate = now.strftime('%Y%m%d')
  rtstruct.StructureSetTime = now.strftime('%H%M%S.%f')
  rtstruct.InstanceCreationDate = now.strftime('%Y%m%d')
  rtstruct.InstanceCreationTime = now.strftime('%H%M%S.%f')
  rtstruct.Manufacturer = '3D Slicer'
  rtstruct.ManufacturerModelName = 'SlicerDynamicPET'
  try:
    import slicer
    rtstruct.SoftwareVersions = str(getattr(slicer.app, 'applicationVersion', ''))
  except Exception:
    pass
  rtstruct.ApprovalStatus = 'UNAPPROVED'

  referenced_frame = Dataset()
  referenced_frame.FrameOfReferenceUID = str(pet_metadata['frameOfReferenceUID'])
  referenced_study = Dataset()
  referenced_study.ReferencedSOPClassUID = RT_REFERENCED_STUDY_SOP_CLASS_UID
  referenced_study.ReferencedSOPInstanceUID = str(pet_metadata['studyInstanceUID'])
  referenced_series = Dataset()
  referenced_series.SeriesInstanceUID = str(pet_metadata['seriesInstanceUID'])
  referenced_series.ContourImageSequence = Sequence([])

  # Complete temporal provenance map. For 3D PET there is one source instance
  # per temporal frame; for classic 2D PET there are multiple slice instances.
  # This top-level map is lightweight and also preserves completely empty
  # temporal segmentation states.
  for frame_index, frame in enumerate(pet_metadata['frames']):
    for instance in frame.get('instances') or []:
      image_reference = Dataset()
      image_reference.ReferencedSOPClassUID = str(instance['sopClassUID'])
      image_reference.ReferencedSOPInstanceUID = str(instance['sopInstanceUID'])
      image_reference.ReferencedFrameNumber = frame_index + 1
      referenced_series.ContourImageSequence.append(image_reference)

  referenced_study.RTReferencedSeriesSequence = Sequence([referenced_series])
  referenced_frame.RTReferencedStudySequence = Sequence([referenced_study])
  rtstruct.ReferencedFrameOfReferenceSequence = Sequence([referenced_frame])
  rtstruct.StructureSetROISequence = Sequence([])
  rtstruct.ROIContourSequence = Sequence([])
  rtstruct.RTROIObservationsSequence = Sequence([])
  return rtstruct

def _frame_plane_candidates(frame):
  cached = frame.get('_planeCandidates')
  if cached is not None:
    return cached
  candidates = []
  for instance in frame.get('instances') or []:
    position = instance.get('imagePositionPatient')
    orientation = instance.get('imageOrientationPatient')
    if not position or not orientation or len(position) != 3 or len(orientation) != 6:
      continue
    row = [float(value) for value in orientation[:3]]
    column = [float(value) for value in orientation[3:]]
    normal = [
      row[1] * column[2] - row[2] * column[1],
      row[2] * column[0] - row[0] * column[2],
      row[0] * column[1] - row[1] * column[0],
    ]
    norm = math.sqrt(sum(value * value for value in normal))
    if norm <= 1e-12:
      continue
    normal = [value / norm for value in normal]
    candidates.append((instance, [float(value) for value in position], normal))
  frame['_planeCandidates'] = candidates
  return candidates


def _reference_instance_for_contour(frame, points_lps):
  """Choose one source PET SOP for a contour without duplicating contours.

    * one-instance/3D frame -> that temporal SOP;
    * multi-instance/2D frame -> nearest source slice plane.
  """
  instances = frame.get('instances') or []
  if len(instances) == 1:
    return instances[0]
  if not instances:
    raise RuntimeError('PET temporal frame has no source SOP instances.')
  candidates = _frame_plane_candidates(frame)
  if not candidates:
    raise RuntimeError(
      'A multi-instance PET temporal frame lacks ImagePositionPatient / '
      'ImageOrientationPatient metadata required to reference its 2D slices.')

  point = points_lps[0]
  best = None
  best_distance = None
  for instance, position, normal in candidates:
    distance = abs(sum(
      (float(point[index]) - position[index]) * normal[index]
      for index in range(3)))
    if best_distance is None or distance < best_distance:
      best = instance
      best_distance = distance
  return best


def _append_roi(rtstruct, roi_number, definition, geometry_states, pet_metadata):
  """Append one semantic ROI using compact unique-geometry serialization.

  geometry_states maps a binary-geometry signature to:
    {'contours': [...], 'frames': [zero-based temporal frame indices]}

  ContourData is written once for each unique geometry.  A contour carries one
  image reference for every temporal state in which that geometry is active.
  This avoids duplicating identical point arrays dozens of times.
  """
  from pydicom.dataset import Dataset
  from pydicom.sequence import Sequence

  structure_set_roi = Dataset()
  structure_set_roi.ROINumber = int(roi_number)
  structure_set_roi.ReferencedFrameOfReferenceUID = str(
    pet_metadata['frameOfReferenceUID'])
  structure_set_roi.ROIName = str(definition['name'])
  structure_set_roi.ROIDescription = 'SlicerDynamicPET dynamic ROI'
  structure_set_roi.ROIGenerationAlgorithm = 'MANUAL'
  rtstruct.StructureSetROISequence.append(structure_set_roi)

  roi_contour = Dataset()
  roi_contour.ReferencedROINumber = int(roi_number)
  roi_contour.ROIDisplayColor = [
    int(round(255.0 * max(0.0, min(1.0, component))))
    for component in definition['color']
  ]
  roi_contour.ContourSequence = Sequence([])

  for state in geometry_states.values():
    contours = state.get('contours') or []
    frame_indices = sorted(set(int(value) for value in (state.get('frames') or [])))
    if not contours or not frame_indices:
      continue

    for points_lps in contours:
      if len(points_lps) < 3:
        continue
      contour = Dataset()
      contour.ContourGeometricType = 'CLOSED_PLANAR'
      contour.NumberOfContourPoints = len(points_lps)
      flattened = []
      for point in points_lps:
        flattened.extend(round(float(value), 5) for value in point)
      contour.ContourData = flattened
      contour.ContourImageSequence = Sequence([])

      seen_references = set()
      for frame_index in frame_indices:
        frame = pet_metadata['frames'][frame_index]
        source_instance = _reference_instance_for_contour(frame, points_lps)
        reference_key = (
          str(source_instance['sopInstanceUID']), frame_index + 1)
        if reference_key in seen_references:
          continue
        seen_references.add(reference_key)
        image_reference = Dataset()
        image_reference.ReferencedSOPClassUID = str(source_instance['sopClassUID'])
        image_reference.ReferencedSOPInstanceUID = str(source_instance['sopInstanceUID'])
        image_reference.ReferencedFrameNumber = frame_index + 1
        contour.ContourImageSequence.append(image_reference)

      if contour.ContourImageSequence:
        roi_contour.ContourSequence.append(contour)

  if not roi_contour.ContourSequence:
    rtstruct.StructureSetROISequence.pop()
    return False

  rtstruct.ROIContourSequence.append(roi_contour)
  observation = Dataset()
  observation.ObservationNumber = int(roi_number)
  observation.ReferencedROINumber = int(roi_number)
  observation.RTROIInterpretedType = ''
  observation.ROIInterpreter = ''
  rtstruct.RTROIObservationsSequence.append(observation)
  return True


def _contour_points_signature(contours):
  """Deterministic digest at the same precision written to DICOM.

  ContourData is serialized after rounding each coordinate to 5 decimal
  places.  Hash the canonical serialized geometry, not the higher-precision
  in-memory coordinates, otherwise valid round-trip geometry can fail the
  identity validator solely because of that intentional quantization.
  """
  digest = hashlib.blake2b(digest_size=16)
  digest.update(struct.pack('<I', len(contours)))
  for points in contours:
    digest.update(struct.pack('<I', len(points)))
    for point in points:
      digest.update(struct.pack(
        '<3d',
        round(float(point[0]), 5),
        round(float(point[1]), 5),
        round(float(point[2]), 5)))
  return digest.hexdigest()


def _dicom_contour_signature_for_frame(roi_contour, frame_number):
  """Digest ContourData belonging to one ROI temporal frame."""
  contours = []
  for contour in getattr(roi_contour, 'ContourSequence', []):
    referenced_frames = {
      int(reference.ReferencedFrameNumber)
      for reference in getattr(contour, 'ContourImageSequence', [])
      if getattr(reference, 'ReferencedFrameNumber', None) is not None
    }
    if int(frame_number) not in referenced_frames:
      continue
    values = [float(value) for value in getattr(contour, 'ContourData', [])]
    points = [
      (values[index], values[index + 1], values[index + 2])
      for index in range(0, len(values), 3)
    ]
    contours.append(points)
  return _contour_points_signature(contours)


def _validate_rtstruct_in_memory(
    rtstruct,
    pet_metadata,
    expected_roi_frame_contour_counts=None,
    expected_roi_frame_geometry_signatures=None):
  """Fast structural validation without re-reading the just-written file."""
  if str(getattr(rtstruct, 'SOPClassUID', '')) != RTSTRUCT_SOP_CLASS_UID:
    raise RuntimeError('Generated object is not an RT Structure Set.')

  expected_map = {
    index + 1: {str(item['sopInstanceUID']) for item in frame.get('instances') or []}
    for index, frame in enumerate(pet_metadata['frames'])
  }
  actual_map = {}
  try:
    series = (rtstruct.ReferencedFrameOfReferenceSequence[0]
              .RTReferencedStudySequence[0]
              .RTReferencedSeriesSequence[0])
    for image_reference in getattr(series, 'ContourImageSequence', []):
      frame_number = int(image_reference.ReferencedFrameNumber)
      actual_map.setdefault(frame_number, set()).add(
        str(image_reference.ReferencedSOPInstanceUID))
  except Exception as error:
    raise RuntimeError('Generated RTSTRUCT is missing the temporal reference map.') from error
  if actual_map != expected_map:
    raise RuntimeError('Generated RTSTRUCT failed temporal SOP-set validation.')

  contour_count = 0
  actual_roi_frame_contour_counts = {}
  for roi_contour in getattr(rtstruct, 'ROIContourSequence', []):
    roi_number = int(getattr(roi_contour, 'ReferencedROINumber', 0) or 0)
    frame_counts = actual_roi_frame_contour_counts.setdefault(roi_number, {})
    for contour in getattr(roi_contour, 'ContourSequence', []):
      contour_count += 1
      if str(getattr(contour, 'ContourGeometricType', '')) != 'CLOSED_PLANAR':
        raise RuntimeError('Generated contour is not CLOSED_PLANAR.')
      point_count = int(getattr(contour, 'NumberOfContourPoints', 0))
      if point_count < 3 or len(getattr(contour, 'ContourData', [])) != 3 * point_count:
        raise RuntimeError('Generated contour has invalid point data.')
      references = getattr(contour, 'ContourImageSequence', [])
      if not references:
        raise RuntimeError('Generated contour has no temporal PET references.')
      referenced_frames_for_contour = set()
      for reference in references:
        frame_number = int(reference.ReferencedFrameNumber)
        uid = str(reference.ReferencedSOPInstanceUID)
        if uid not in expected_map.get(frame_number, set()):
          raise RuntimeError('Generated contour failed temporal PET reference validation.')
        referenced_frames_for_contour.add(frame_number)
      for frame_number in referenced_frames_for_contour:
        frame_counts[frame_number] = frame_counts.get(frame_number, 0) + 1
  if contour_count == 0:
    raise RuntimeError('Generated RTSTRUCT contains no ROI contours.')
  if expected_roi_frame_contour_counts is not None:
    normalized_actual = {
      int(roi): {int(frame): int(count) for frame, count in counts.items()}
      for roi, counts in actual_roi_frame_contour_counts.items()
    }
    normalized_expected = {
      int(roi): {int(frame): int(count) for frame, count in counts.items()}
      for roi, counts in expected_roi_frame_contour_counts.items()
    }
    if normalized_actual != normalized_expected:
      raise RuntimeError(
        'Generated RTSTRUCT failed per-ROI temporal contour-state validation.')

  if expected_roi_frame_geometry_signatures is not None:
    actual_geometry = {}
    for roi_contour in getattr(rtstruct, 'ROIContourSequence', []):
      roi_number = int(getattr(roi_contour, 'ReferencedROINumber', 0) or 0)
      frame_numbers = sorted(
        actual_roi_frame_contour_counts.get(roi_number, {}).keys())
      actual_geometry[roi_number] = {
        int(frame_number): _dicom_contour_signature_for_frame(
          roi_contour, frame_number)
        for frame_number in frame_numbers
      }
    normalized_expected_geometry = {
      int(roi): {int(frame): str(signature) for frame, signature in states.items()}
      for roi, states in expected_roi_frame_geometry_signatures.items()
    }
    if actual_geometry != normalized_expected_geometry:
      mismatch_details = []
      roi_numbers = sorted(
        set(actual_geometry.keys()) | set(normalized_expected_geometry.keys()))
      for roi_number in roi_numbers:
        actual_frames = actual_geometry.get(roi_number, {})
        expected_frames = normalized_expected_geometry.get(roi_number, {})
        frame_numbers = sorted(set(actual_frames.keys()) | set(expected_frames.keys()))
        for frame_number in frame_numbers:
          actual_signature = actual_frames.get(frame_number)
          expected_signature = expected_frames.get(frame_number)
          if actual_signature != expected_signature:
            mismatch_details.append(
              f'ROI {roi_number} frame {frame_number}: '
              f'expected={expected_signature} actual={actual_signature}')
            if len(mismatch_details) >= 5:
              break
        if len(mismatch_details) >= 5:
          break
      detail_text = '; '.join(mismatch_details) or 'unknown mismatch'
      raise RuntimeError(
        'Generated RTSTRUCT failed per-ROI temporal geometry identity validation. '
        + detail_text)

  return contour_count

def _format_frame_ranges(zero_based_frames):
  values = sorted(set(int(value) + 1 for value in zero_based_frames))
  if not values:
    return '-'
  ranges = []
  start = previous = values[0]
  for value in values[1:]:
    if value == previous + 1:
      previous = value
      continue
    ranges.append(str(start) if start == previous else f'{start}-{previous}')
    start = previous = value
  ranges.append(str(start) if start == previous else f'{start}-{previous}')
  return ','.join(ranges)


def export_dynamic_rtstruct(
    segmentation_sequence_node,
    pet_sequence_node,
    reference_volume_node,
    output_path,
    *,
    pet_dicom_files=None,
    series_description=None,
    series_number=301,
    structure_set_label='DynamicRT',
    require_enhanced_pet=False,
    overwrite=False,
    progress_callback=None):
  """Export a compact temporal RTSTRUCT.

  PET provenance priority:
    1. dPETImporter persisted frame references;
    2. MRML + Slicer DICOM database metadata;
    3. source DICOM headers as a final fallback.

  Identical segment geometries are contoured and serialized once. Their
  ContourImageSequence contains all temporal PET references for which that
  geometry applies.
  """
  import pydicom

  start_time = time.perf_counter()
  if segmentation_sequence_node is None:
    raise ValueError('segmentation_sequence_node is required.')
  if pet_sequence_node is None:
    raise ValueError('pet_sequence_node is required.')
  if segmentation_sequence_node.GetNumberOfDataNodes() < 1:
    raise ValueError('The segmentation sequence is empty.')
  if pet_sequence_node.GetNumberOfDataNodes() < 1:
    raise ValueError('The PET sequence is empty.')

  pet_frame_count = pet_sequence_node.GetNumberOfDataNodes()
  if segmentation_sequence_node.GetNumberOfDataNodes() > pet_frame_count:
    raise ValueError('The segmentation sequence contains more frames than the PET sequence.')

  output_path = os.path.abspath(str(output_path))
  if not output_path.lower().endswith('.dcm'):
    output_path += '.dcm'
  if os.path.exists(output_path) and not overwrite:
    raise FileExistsError(f'Output already exists: {output_path}')

  def report(percent, text):
    if progress_callback is None:
      return
    if progress_callback(int(percent), 100, text) is False:
      raise RuntimeError('Dynamic RTSTRUCT export was cancelled.')

  report(0, 'Resolving PET DICOM metadata')
  pet_metadata = _resolve_pet_reference_metadata(
    pet_sequence_node, pet_dicom_files, require_enhanced_pet)

  reference_frame_uid = _reference_frame_of_reference_uid(reference_volume_node)
  pet_frame_uid = str(pet_metadata['frameOfReferenceUID'])
  if reference_frame_uid and reference_frame_uid != pet_frame_uid:
    raise ValueError(
      'The selected export geometry reference and dynamic PET have different '
      'DICOM Frame of Reference UIDs.')

  report(5, 'Analyzing segmentation')
  definitions = _segment_definitions(segmentation_sequence_node)
  if not definitions:
    raise ValueError('No segments are present in the segmentation sequence.')

  # Per semantic segment, keep only unique geometries and the frames in which
  # each geometry is active.  Keep the contour cache ROI-scoped as well: even
  # if two segments happen to have byte-identical masks, no mutable contour
  # container is shared across ROI identities.
  temporal_geometries = {segment_id: {} for segment_id in definitions}
  contour_cache = {}
  cache_hits = 0
  cache_misses = 0
  extraction_start = time.perf_counter()

  for frame_index in range(pet_frame_count):
    index_value = pet_sequence_node.GetNthIndexValue(frame_index)
    frame_segmentation = segmentation_sequence_node.GetDataNodeAtValue(
      index_value, True)

    if frame_segmentation is not None:
      frame_segment_ids = set(_segment_ids(frame_segmentation))
      for segment_id in definitions:
        if segment_id not in frame_segment_ids:
          continue
        labelmap = _binary_labelmap(frame_segmentation, segment_id)
        if labelmap is None:
          continue

        signature = _labelmap_signature(labelmap)
        cache_key = (segment_id, signature)
        if cache_key in contour_cache:
          contours = contour_cache[cache_key]
          cache_hits += 1
        else:
          contours = _contours_from_binary_labelmap(labelmap)
          contour_cache[cache_key] = contours
          cache_misses += 1

        if not contours:
          continue
        state = temporal_geometries[segment_id].get(signature)
        if state is None:
          state = {'contours': contours, 'frames': []}
          temporal_geometries[segment_id][signature] = state
        state['frames'].append(frame_index)

    percent = 5 + int(round(55.0 * (frame_index + 1) / pet_frame_count))
    report(percent, f'Extracting unique contours: frame {frame_index + 1}/{pet_frame_count}')

  # Explicitly publish the terminal frame count before moving to the next
  # phase.  On fast final frames Qt could repaint directly from 65/66 to the
  # assembly phase, making the extraction counter look off by one.
  report(60, f'Contours extracted: {pet_frame_count}/{pet_frame_count} frames')

  extraction_seconds = time.perf_counter() - extraction_start
  unique_states = sum(len(states) for states in temporal_geometries.values())
  logging.info(
    '[dRTExporter] Contour extraction %.3f s; cache hits=%d misses=%d; unique ROI states=%d',
    extraction_seconds, cache_hits, cache_misses, unique_states)
  for segment_id, definition in definitions.items():
    states = temporal_geometries.get(segment_id, {})
    state_text = '; '.join(
      f"state={signature[:10]} {len(state.get('contours') or [])} contours @ frames "
      f"{_format_frame_ranges(state.get('frames') or [])}"
      for signature, state in states.items()) or 'empty'
    logging.info(
      "[dRTExporter] ROI '%s' id=%s: %d unique state(s): %s",
      definition['name'], segment_id, len(states), state_text)

  report(65, 'Building compact dynamic RTSTRUCT')
  assembly_start = time.perf_counter()

  if not series_description:
    series_description = (
      segmentation_sequence_node.GetAttribute('dRTImporter.DisplayName')
      or segmentation_sequence_node.GetName()
      or 'SlicerDynamicPET')
    # Keep the DICOM series name concise when exporting a sequence node whose
    # MRML name already carries an implementation suffix.
    for suffix in (' [dRT]', ' segmentation sequence', ' sequence'):
      if series_description.lower().endswith(suffix.lower()):
        series_description = series_description[:-len(suffix)].strip()
        break
    if not series_description.lower().endswith('dynamic rtstruct'):
      series_description = f'{series_description} - Dynamic RTSTRUCT'

  rtstruct = _new_rtstruct_dataset(
    pet_metadata,
    series_description=series_description,
    series_number=series_number,
    structure_set_label=structure_set_label)

  exported_roi_number = 0
  expected_roi_frame_contour_counts = {}
  expected_roi_frame_geometry_signatures = {}
  for segment_id, definition in definitions.items():
    geometry_states = temporal_geometries.get(segment_id, {})
    if not geometry_states:
      continue
    proposed_roi_number = exported_roi_number + 1
    if _append_roi(
        rtstruct, proposed_roi_number, definition, geometry_states, pet_metadata):
      exported_roi_number = proposed_roi_number
      frame_counts = {}
      for state in geometry_states.values():
        contour_count_for_state = len(state.get('contours') or [])
        for frame_index in state.get('frames') or []:
          frame_number = int(frame_index) + 1
          if frame_number in frame_counts:
            raise RuntimeError(
              f"ROI '{definition['name']}' was assigned more than one geometry "
              f'to temporal frame {frame_number}.')
          frame_counts[frame_number] = contour_count_for_state
      expected_roi_frame_contour_counts[proposed_roi_number] = frame_counts
      frame_geometry_signatures = {}
      for state in geometry_states.values():
        geometry_signature = _contour_points_signature(state.get('contours') or [])
        for frame_index in state.get('frames') or []:
          frame_number = int(frame_index) + 1
          frame_geometry_signatures[frame_number] = geometry_signature
      expected_roi_frame_geometry_signatures[proposed_roi_number] = (
        frame_geometry_signatures)

  if exported_roi_number == 0:
    raise ValueError('No non-empty ROI contours were generated.')

  contour_count = _validate_rtstruct_in_memory(
    rtstruct,
    pet_metadata,
    expected_roi_frame_contour_counts,
    expected_roi_frame_geometry_signatures)
  assembly_seconds = time.perf_counter() - assembly_start
  logging.info(
    '[dRTExporter] Compact RTSTRUCT assembly %.3f s; ROIs=%d contours=%d',
    assembly_seconds, exported_roi_number, contour_count)

  report(85, 'Writing dynamic RTSTRUCT DICOM')
  os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
  write_start = time.perf_counter()
  pydicom.dcmwrite(output_path, rtstruct, enforce_file_format=True)
  write_seconds = time.perf_counter() - write_start
  if not os.path.isfile(output_path) or os.path.getsize(output_path) <= 132:
    raise RuntimeError('Dynamic RTSTRUCT write did not create a valid DICOM file.')

  report(100, 'Dynamic RTSTRUCT export complete')
  total_seconds = time.perf_counter() - start_time
  size_mb = os.path.getsize(output_path) / (1024.0 * 1024.0)
  logging.info(
    '[dRTExporter] Export complete %.3f s; write %.3f s; size %.2f MB: %s',
    total_seconds, write_seconds, size_mb, output_path)
  return output_path

def export_dynamic_rtstruct_from_node_ids(
    segmentation_sequence_node_id,
    pet_sequence_node_id,
    reference_volume_node_id,
    output_path,
    overwrite=False,
    show_progress=True):
  """PythonQt-friendly wrapper used by the DynamicPET C++ widget."""
  import qt
  import slicer

  segmentation_sequence_node = slicer.mrmlScene.GetNodeByID(
    str(segmentation_sequence_node_id))
  pet_sequence_node = slicer.mrmlScene.GetNodeByID(str(pet_sequence_node_id))
  reference_volume_node = (
    slicer.mrmlScene.GetNodeByID(str(reference_volume_node_id))
    if reference_volume_node_id else None)

  total_steps = 100

  progress = None
  if show_progress:
    progress = slicer.util.createProgressDialog(
      labelText='Preparing dynamic RTSTRUCT export',
      value=0,
      maximum=total_steps,
      windowModality=qt.Qt.WindowModal)

  def progress_callback(value, maximum, label):
    if progress is None:
      return True
    progress.maximum = maximum
    progress.value = value
    progress.labelText = label
    slicer.app.processEvents()
    return not progress.wasCanceled

  _start_diagnostics()
  try:
    _dbg(
      f'wrapper nodeIDs segSeq={segmentation_sequence_node_id!r} '
      f'petSeq={pet_sequence_node_id!r} ref={reference_volume_node_id!r}')
    return export_dynamic_rtstruct(
      segmentation_sequence_node,
      pet_sequence_node,
      reference_volume_node,
      output_path,
      overwrite=bool(overwrite),
      progress_callback=progress_callback)
  finally:
    if progress is not None:
      progress.close()
    _stop_diagnostics()


__all__ = [
  'export_dynamic_rtstruct',
  'export_dynamic_rtstruct_from_node_ids',
]
