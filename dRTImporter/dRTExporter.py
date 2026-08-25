"""Dynamic RT Structure Set exporter for SlicerDynamicPET.

The exporter runs inside 3D Slicer. It uses SlicerRT to convert each temporal
segmentation state into a conventional RT Structure Set, then uses pydicom to
assemble those contours into one temporal RTSTRUCT convention understood by
``dRTImporter``.

No contour approximation package is used: contour generation is delegated to
SlicerRT, which is already required by the importer.
"""

import copy
import hashlib
import json
import os
import tempfile
from datetime import datetime


ENHANCED_PET_IMAGE_STORAGE_UID = "1.2.840.10008.5.1.4.1.1.130"
RTSTRUCT_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.481.3"


def _segment_ids(segmentation_node):
  import vtk

  ids = vtk.vtkStringArray()
  segmentation_node.GetSegmentation().GetSegmentIDs(ids)
  return [ids.GetValue(index) for index in range(ids.GetNumberOfValues())]


def _pet_source_files(pet_sequence_node, explicit_files=None):
  """Resolve exactly one source DICOM object for each PET temporal frame."""
  if explicit_files is not None:
    files = [os.path.abspath(str(path)) for path in explicit_files]
  else:
    frame_file_list = pet_sequence_node.GetAttribute('MultiVolume.FrameFileList')
    if frame_file_list:
      files = [path for path in frame_file_list.split(',') if path]
    else:
      files = []
      for frame_index in range(pet_sequence_node.GetNumberOfDataNodes()):
        frame_node = pet_sequence_node.GetNthDataNode(frame_index)
        source_files_json = frame_node.GetAttribute('dPET.SourceDICOMFiles')
        if not source_files_json:
          raise ValueError(
            f'PET frame {frame_index + 1} has no source DICOM file metadata.')
        frame_files = json.loads(source_files_json)
        if len(frame_files) != 1:
          raise ValueError(
            'Dynamic RTSTRUCT export currently requires one source DICOM '
            f'object per temporal PET frame; frame {frame_index + 1} has '
            f'{len(frame_files)} source files.')
        files.append(frame_files[0])

  expected_count = pet_sequence_node.GetNumberOfDataNodes()
  if len(files) != expected_count:
    raise ValueError(
      'Dynamic RTSTRUCT export requires one source DICOM object per temporal '
      f'PET frame; found {len(files)} files for {expected_count} frames.')

  missing = [path for path in files if not os.path.isfile(path)]
  if missing:
    raise FileNotFoundError(f'PET source DICOM file is unavailable: {missing[0]}')
  return files


def _validate_pet_datasets(pet_datasets, expected_count, require_enhanced_pet):
  if len(pet_datasets) != expected_count:
    raise ValueError('PET dataset count does not match the PET sequence.')

  sop_instance_uids = [str(dataset.SOPInstanceUID) for dataset in pet_datasets]
  if len(set(sop_instance_uids)) != len(sop_instance_uids):
    raise ValueError(
      'Each temporal PET frame must reference a distinct SOP Instance UID.')

  study_uids = {str(getattr(dataset, 'StudyInstanceUID', '')) for dataset in pet_datasets}
  series_uids = {str(getattr(dataset, 'SeriesInstanceUID', '')) for dataset in pet_datasets}
  frame_of_reference_uids = {
    str(getattr(dataset, 'FrameOfReferenceUID', '')) for dataset in pet_datasets
  }
  if len(study_uids) != 1 or '' in study_uids:
    raise ValueError('PET frames must have one valid Study Instance UID.')
  if len(series_uids) != 1 or '' in series_uids:
    raise ValueError('PET frames must have one valid Series Instance UID.')
  if len(frame_of_reference_uids) != 1 or '' in frame_of_reference_uids:
    raise ValueError('PET frames must have one valid Frame of Reference UID.')

  if require_enhanced_pet:
    non_enhanced = [
      str(dataset.SOPClassUID)
      for dataset in pet_datasets
      if str(dataset.SOPClassUID) != ENHANCED_PET_IMAGE_STORAGE_UID
    ]
    if non_enhanced:
      raise ValueError(
        'Dynamic RTSTRUCT export expects Enhanced PET Image Storage when '
        f'require_enhanced_pet=True. Found SOP Class UID {non_enhanced[0]}.')


def _reference_frame_of_reference_uid(reference_volume_node):
  """Resolve the source DICOM Frame of Reference UID when available."""
  import slicer
  import pydicom

  direct = reference_volume_node.GetAttribute('DICOM.FrameOfReferenceUID')
  if direct:
    return str(direct)

  instance_uids = (reference_volume_node.GetAttribute('DICOM.instanceUIDs') or '').split()
  if not instance_uids:
    return None
  try:
    file_path = slicer.dicomDatabase.fileForInstance(instance_uids[0])
    if not file_path or not os.path.isfile(file_path):
      return None
    dataset = pydicom.dcmread(
      file_path, stop_before_pixels=True, force=True,
      specific_tags=['FrameOfReferenceUID'])
    value = str(getattr(dataset, 'FrameOfReferenceUID', '') or '')
    return value or None
  except Exception:
    return None


def _stable_roi_token(segment_id):
  digest = hashlib.sha1(segment_id.encode('utf-8')).hexdigest()[:16]
  return f'SDP_{digest}'


def _segment_definitions(segmentation_sequence_node):
  """Return stable segment definitions keyed by persistent Slicer segment ID."""
  definitions = {}
  for item_index in range(segmentation_sequence_node.GetNumberOfDataNodes()):
    segmentation_node = segmentation_sequence_node.GetNthDataNode(item_index)
    for segment_id in _segment_ids(segmentation_node):
      segment = segmentation_node.GetSegmentation().GetSegment(segment_id)
      if segment_id not in definitions:
        definitions[segment_id] = {
          'id': segment_id,
          'name': segment.GetName(),
          'color': tuple(float(value) for value in segment.GetColor()),
          'token': _stable_roi_token(segment_id),
        }
  return definitions


def _nonempty_segment_ids(segmentation_node, reference_volume_node):
  """Return segment IDs that contain at least one voxel in reference geometry."""
  import numpy as np
  import slicer

  nonempty = []
  for segment_id in _segment_ids(segmentation_node):
    try:
      labelmap = slicer.util.arrayFromSegmentBinaryLabelmap(
        segmentation_node, segment_id, reference_volume_node)
    except Exception:
      # Let SlicerRT attempt conversion if the binary representation is not
      # immediately available in the requested geometry.
      nonempty.append(segment_id)
      continue
    if labelmap.size and np.any(labelmap):
      nonempty.append(segment_id)
  return nonempty


def _temporary_segmentation_for_all_frames(
    segmentation_sequence_node,
    pet_sequence_node,
    segment_definitions,
    reference_volume_node,
    progress_callback=None):
  """Create one temporary segmentation containing frame-tagged ROIs.

  SlicerRT then needs to export the reference image series only once. Each
  non-empty (segment, temporal frame) pair becomes one temporary ROI with a
  short unique token. The caller later regroups the resulting contours by the
  persistent Slicer segment ID.
  """
  import slicer
  import vtkSegmentationCorePython as vtkSegmentationCore

  temporary_node = slicer.mrmlScene.AddNewNodeByClass(
    'vtkMRMLSegmentationNode', 'SlicerDynamicPET temporary dynamic RTSTRUCT')
  temporary_node.CreateDefaultDisplayNodes()
  temporary_node.SetReferenceImageGeometryParameterFromVolumeNode(reference_volume_node)

  target_segmentation = temporary_node.GetSegmentation()
  token_map = {}
  pet_frame_count = pet_sequence_node.GetNumberOfDataNodes()

  for pet_frame_index in range(pet_frame_count):
    if progress_callback is not None:
      if progress_callback(pet_frame_index, pet_frame_count) is False:
        slicer.mrmlScene.RemoveNode(temporary_node)
        raise RuntimeError('Dynamic RTSTRUCT export was cancelled.')

    index_value = pet_sequence_node.GetNthIndexValue(pet_frame_index)
    frame_segmentation = segmentation_sequence_node.GetDataNodeAtValue(
      index_value, True)
    if frame_segmentation is None:
      continue

    nonempty_ids = _nonempty_segment_ids(
      frame_segmentation, reference_volume_node)
    source_segmentation = frame_segmentation.GetSegmentation()
    for segment_id in nonempty_ids:
      source_segment = source_segmentation.GetSegment(segment_id)
      definition = segment_definitions.get(segment_id)
      if source_segment is None or definition is None:
        continue

      token = f"{definition['token']}_F{pet_frame_index + 1:06d}"
      segment_copy = vtkSegmentationCore.vtkSegment()
      segment_copy.DeepCopy(source_segment)
      segment_copy.SetName(token)
      temporary_id = f"dRT_{definition['token']}_{pet_frame_index + 1:06d}"
      target_segmentation.AddSegment(segment_copy, temporary_id)
      token_map[token] = {
        'segment_id': segment_id,
        'pet_frame_index': pet_frame_index,
      }

  return temporary_node, token_map

def _slicerrt_export_static_rtstruct(segmentation_node, reference_volume_node, output_directory):
  """Export one conventional RTSTRUCT using the installed SlicerRT plugin."""
  import pydicom
  import slicer
  from DicomRtImportExportPlugin import DicomRtImportExportPluginClass

  if not hasattr(slicer.modules, 'dicomrtimportexport'):
    raise RuntimeError(
      'SlicerRT is required. Install the SlicerRT extension and restart 3D Slicer.')

  sh_node = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
  if sh_node is None:
    raise RuntimeError('Subject hierarchy is unavailable.')

  reference_item_id = sh_node.GetItemByDataNode(reference_volume_node)
  segmentation_item_id = sh_node.GetItemByDataNode(segmentation_node)
  if reference_item_id == slicer.vtkMRMLSubjectHierarchyNode.INVALID_ITEM_ID:
    raise RuntimeError('The reference volume is not present in subject hierarchy.')
  if segmentation_item_id == slicer.vtkMRMLSubjectHierarchyNode.INVALID_ITEM_ID:
    raise RuntimeError('The temporary segmentation is not present in subject hierarchy.')

  reference_parent_id = sh_node.GetItemParent(reference_item_id)
  if reference_parent_id != slicer.vtkMRMLSubjectHierarchyNode.INVALID_ITEM_ID:
    sh_node.SetItemParent(segmentation_item_id, reference_parent_id)

  plugin = DicomRtImportExportPluginClass()
  reference_exportables = plugin.examineForExport(reference_item_id)
  segmentation_exportables = plugin.examineForExport(segmentation_item_id)
  if not reference_exportables:
    raise RuntimeError('SlicerRT could not create a DICOM exportable for the reference volume.')
  if not segmentation_exportables:
    raise RuntimeError('SlicerRT could not create a DICOM RTSTRUCT exportable for the segmentation.')

  reference_exportable = reference_exportables[0]
  segmentation_exportable = segmentation_exportables[0]
  for exportable in (reference_exportable, segmentation_exportable):
    exportable.directory = output_directory
  reference_exportable.setTag('SeriesDescription', 'SlicerDynamicPET temporary reference')
  reference_exportable.setTag('SeriesNumber', '990')
  segmentation_exportable.setTag('SeriesDescription', 'SlicerDynamicPET temporary RTSTRUCT')
  segmentation_exportable.setTag('SeriesNumber', '991')

  message = plugin.export([reference_exportable, segmentation_exportable])
  if message:
    raise RuntimeError(f'SlicerRT RTSTRUCT export failed: {message}')

  rtstruct_paths = []
  for root, _dirs, files in os.walk(output_directory):
    for file_name in files:
      path = os.path.join(root, file_name)
      try:
        dataset = pydicom.dcmread(path, stop_before_pixels=True, force=True)
      except Exception:
        continue
      if (str(getattr(dataset, 'SOPClassUID', '')) == RTSTRUCT_SOP_CLASS_UID
          or str(getattr(dataset, 'Modality', '')) == 'RTSTRUCT'):
        rtstruct_paths.append(path)

  if len(rtstruct_paths) != 1:
    raise RuntimeError(
      'SlicerRT temporary export did not produce exactly one RTSTRUCT '
      f'(found {len(rtstruct_paths)}).')
  return rtstruct_paths[0]


def _roi_contours_by_token(static_rtstruct):
  """Return token -> copied contour datasets from a SlicerRT RTSTRUCT."""
  roi_name_by_number = {
    int(item.ROINumber): str(item.ROIName)
    for item in getattr(static_rtstruct, 'StructureSetROISequence', [])
  }
  result = {}
  for roi_contour in getattr(static_rtstruct, 'ROIContourSequence', []):
    roi_number = int(roi_contour.ReferencedROINumber)
    token = roi_name_by_number.get(roi_number)
    if not token:
      continue
    result[token] = [copy.deepcopy(contour) for contour in getattr(roi_contour, 'ContourSequence', [])]
  return result


def _build_pet_references(rtstruct, pet_datasets):
  """Replace the temporary reference series with the dynamic PET evidence."""
  from pydicom.dataset import Dataset
  from pydicom.sequence import Sequence

  pet0 = pet_datasets[0]
  frame_of_reference_uid = str(pet0.FrameOfReferenceUID)

  referenced_frame = Dataset()
  referenced_frame.FrameOfReferenceUID = frame_of_reference_uid

  referenced_study = Dataset()
  referenced_study.ReferencedSOPClassUID = '1.2.840.10008.3.1.2.3.2'
  referenced_study.ReferencedSOPInstanceUID = str(pet0.StudyInstanceUID)

  referenced_series = Dataset()
  referenced_series.SeriesInstanceUID = str(pet0.SeriesInstanceUID)
  referenced_series.ContourImageSequence = Sequence([])

  # Store the complete temporal map at series level. This lets dRTImporter
  # preserve intentionally empty temporal segmentation states as well.
  for frame_index, pet_dataset in enumerate(pet_datasets):
    image_reference = Dataset()
    image_reference.ReferencedSOPClassUID = str(pet_dataset.SOPClassUID)
    image_reference.ReferencedSOPInstanceUID = str(pet_dataset.SOPInstanceUID)
    image_reference.ReferencedFrameNumber = frame_index + 1
    referenced_series.ContourImageSequence.append(image_reference)

  referenced_study.RTReferencedSeriesSequence = Sequence([referenced_series])
  referenced_frame.RTReferencedStudySequence = Sequence([referenced_study])
  rtstruct.ReferencedFrameOfReferenceSequence = Sequence([referenced_frame])


def export_dynamic_rtstruct(
    segmentation_sequence_node,
    pet_sequence_node,
    reference_volume_node,
    output_path,
    *,
    pet_dicom_files=None,
    series_description='SlicerDynamicPET dynamic RTSTRUCT',
    series_number=301,
    structure_set_label='DynamicRT',
    require_enhanced_pet=True,
    overwrite=False,
    progress_callback=None):
  """Export one segmentation sequence as a dynamic RT Structure Set.

  The temporal convention is intentionally explicit: every contour references
  the PET SOP Instance for its temporal position and carries a one-based
  Referenced Frame Number. A complete frame-to-SOP map is also written in the
  top-level referenced series so empty temporal segmentation states survive a
  SlicerDynamicPET export/import round trip.
  """
  import pydicom
  import slicer
  from pydicom.dataset import Dataset
  from pydicom.sequence import Sequence
  from pydicom.uid import generate_uid

  if segmentation_sequence_node is None:
    raise ValueError('segmentation_sequence_node is required.')
  if pet_sequence_node is None:
    raise ValueError('pet_sequence_node is required.')
  if reference_volume_node is None:
    raise ValueError('reference_volume_node is required.')
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

  pet_files = _pet_source_files(pet_sequence_node, explicit_files=pet_dicom_files)
  pet_datasets = [
    pydicom.dcmread(path, stop_before_pixels=True, force=True) for path in pet_files
  ]
  _validate_pet_datasets(pet_datasets, pet_frame_count, require_enhanced_pet)

  reference_frame_uid = _reference_frame_of_reference_uid(reference_volume_node)
  pet_frame_uid = str(pet_datasets[0].FrameOfReferenceUID)
  if reference_frame_uid and reference_frame_uid != pet_frame_uid:
    raise ValueError(
      'The selected export geometry reference and dynamic PET have different '
      'DICOM Frame of Reference UIDs.')

  definitions = _segment_definitions(segmentation_sequence_node)
  if not definitions:
    raise ValueError('No segments are present in the segmentation sequence.')

  contours_by_segment = {segment_id: [] for segment_id in definitions}

  with tempfile.TemporaryDirectory(prefix='SlicerDynamicPET_dRT_export_') as temporary_root:
    temporary_segmentation, token_map = _temporary_segmentation_for_all_frames(
      segmentation_sequence_node,
      pet_sequence_node,
      definitions,
      reference_volume_node,
      progress_callback=progress_callback)
    try:
      if not token_map:
        raise ValueError(
          'No non-empty contours were generated from the segmentation sequence.')
      static_path = _slicerrt_export_static_rtstruct(
        temporary_segmentation, reference_volume_node, temporary_root)
      static_rtstruct = pydicom.dcmread(static_path, force=True)
    finally:
      slicer.mrmlScene.RemoveNode(temporary_segmentation)

  rtstruct = copy.deepcopy(static_rtstruct)
  frame_contours = _roi_contours_by_token(static_rtstruct)
  for token, contours in frame_contours.items():
    mapping = token_map.get(token)
    if mapping is None:
      continue
    segment_id = mapping['segment_id']
    pet_frame_index = mapping['pet_frame_index']
    pet_dataset = pet_datasets[pet_frame_index]
    for contour in contours:
      image_reference = Dataset()
      image_reference.ReferencedSOPClassUID = str(pet_dataset.SOPClassUID)
      image_reference.ReferencedSOPInstanceUID = str(pet_dataset.SOPInstanceUID)
      image_reference.ReferencedFrameNumber = pet_frame_index + 1
      contour.ContourImageSequence = Sequence([image_reference])
      contours_by_segment[segment_id].append(contour)

  rtstruct.StructureSetROISequence = Sequence([])
  rtstruct.ROIContourSequence = Sequence([])
  rtstruct.RTROIObservationsSequence = Sequence([])
  _build_pet_references(rtstruct, pet_datasets)

  frame_of_reference_uid = str(pet_datasets[0].FrameOfReferenceUID)
  exported_roi_number = 0
  for segment_id, definition in definitions.items():
    contours = contours_by_segment.get(segment_id, [])
    if not contours:
      continue

    exported_roi_number += 1
    structure_set_roi = Dataset()
    structure_set_roi.ROINumber = exported_roi_number
    structure_set_roi.ReferencedFrameOfReferenceUID = frame_of_reference_uid
    structure_set_roi.ROIName = str(definition['name'])
    structure_set_roi.ROIDescription = 'SlicerDynamicPET dynamic ROI'
    structure_set_roi.ROIGenerationAlgorithm = 'MANUAL'
    rtstruct.StructureSetROISequence.append(structure_set_roi)

    roi_contour = Dataset()
    roi_contour.ReferencedROINumber = exported_roi_number
    roi_contour.ROIDisplayColor = [
      int(round(255.0 * max(0.0, min(1.0, component))))
      for component in definition['color']
    ]
    roi_contour.ContourSequence = Sequence(contours)
    rtstruct.ROIContourSequence.append(roi_contour)

    observation = Dataset()
    observation.ObservationNumber = exported_roi_number
    observation.ReferencedROINumber = exported_roi_number
    observation.RTROIInterpretedType = ''
    observation.ROIInterpreter = ''
    rtstruct.RTROIObservationsSequence.append(observation)

  if not rtstruct.ROIContourSequence:
    raise ValueError('No non-empty ROI contours were generated.')

  now = datetime.now()
  sop_instance_uid = generate_uid()
  rtstruct.SOPInstanceUID = sop_instance_uid
  if getattr(rtstruct, 'file_meta', None) is not None:
    rtstruct.file_meta.MediaStorageSOPInstanceUID = sop_instance_uid
  rtstruct.SeriesInstanceUID = generate_uid()
  rtstruct.StudyInstanceUID = str(pet_datasets[0].StudyInstanceUID)
  rtstruct.SeriesDescription = str(series_description)
  rtstruct.SeriesNumber = int(series_number)
  rtstruct.StructureSetLabel = str(structure_set_label)[:16]
  rtstruct.StructureSetName = str(series_description)
  rtstruct.StructureSetDescription = 'Dynamic RTSTRUCT generated by 3D Slicer / SlicerDynamicPET'
  rtstruct.StructureSetDate = now.strftime('%Y%m%d')
  rtstruct.StructureSetTime = now.strftime('%H%M%S.%f')
  rtstruct.SeriesDate = now.strftime('%Y%m%d')
  rtstruct.SeriesTime = now.strftime('%H%M%S.%f')
  rtstruct.InstanceCreationDate = now.strftime('%Y%m%d')
  rtstruct.InstanceCreationTime = now.strftime('%H%M%S.%f')
  rtstruct.Manufacturer = '3D Slicer'
  rtstruct.ManufacturerModelName = 'SlicerDynamicPET'
  rtstruct.SoftwareVersions = str(getattr(slicer.app, 'applicationVersion', ''))
  rtstruct.ApprovalStatus = 'UNAPPROVED'

  os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
  pydicom.dcmwrite(output_path, rtstruct, write_like_original=False)

  # Structural round-trip validation.
  saved = pydicom.dcmread(output_path, stop_before_pixels=True, force=True)
  top_level_map = {}
  try:
    series = (saved.ReferencedFrameOfReferenceSequence[0]
              .RTReferencedStudySequence[0]
              .RTReferencedSeriesSequence[0])
    for image_reference in getattr(series, 'ContourImageSequence', []):
      frame_number = int(image_reference.ReferencedFrameNumber)
      top_level_map[frame_number] = str(image_reference.ReferencedSOPInstanceUID)
  except Exception as error:
    raise RuntimeError('Saved RTSTRUCT is missing the temporal reference map.') from error

  expected_map = {
    index + 1: str(dataset.SOPInstanceUID)
    for index, dataset in enumerate(pet_datasets)
  }
  if top_level_map != expected_map:
    raise RuntimeError('Saved RTSTRUCT failed temporal SOP reference validation.')

  for roi_contour in saved.ROIContourSequence:
    for contour in roi_contour.ContourSequence:
      references = getattr(contour, 'ContourImageSequence', [])
      if len(references) != 1:
        raise RuntimeError('Saved contour does not have exactly one temporal PET reference.')
      reference = references[0]
      frame_number = int(reference.ReferencedFrameNumber)
      if expected_map.get(frame_number) != str(reference.ReferencedSOPInstanceUID):
        raise RuntimeError('Saved contour failed temporal PET reference validation.')

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

  segmentation_sequence_node = slicer.mrmlScene.GetNodeByID(str(segmentation_sequence_node_id))
  pet_sequence_node = slicer.mrmlScene.GetNodeByID(str(pet_sequence_node_id))
  reference_volume_node = slicer.mrmlScene.GetNodeByID(str(reference_volume_node_id))

  progress = None
  if show_progress:
    progress = slicer.util.createProgressDialog(
      labelText='Exporting dynamic RTSTRUCT',
      value=0,
      maximum=max(1, pet_sequence_node.GetNumberOfDataNodes() if pet_sequence_node else 1),
      windowModality=qt.Qt.WindowModal)

  def progress_callback(frame_index, frame_count):
    if progress is None:
      return True
    progress.value = frame_index
    progress.maximum = frame_count
    progress.labelText = f'Exporting dynamic RTSTRUCT frame {frame_index + 1}/{frame_count}'
    slicer.app.processEvents()
    return not progress.wasCanceled

  try:
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


__all__ = [
  'export_dynamic_rtstruct',
  'export_dynamic_rtstruct_from_node_ids',
]
