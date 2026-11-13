"""Strict support for SlicerDynamicPET's Enhanced single-SOP dynamic RTSTRUCT.

This convention is specific to DynamicPET.  It is intentionally not inferred
from ordinary RTSTRUCT ReferencedFrameNumber, which denotes an actual DICOM
image frame, not a reconstructed PET temporal volume.
"""

import json
from numbers import Integral

ENHANCED_PET_SOP_CLASS_UID = '1.2.840.10008.5.1.4.1.1.130'
RTSTRUCT_SOP_CLASS_UID = '1.2.840.10008.5.1.4.1.1.481.3'


def _numbers(value):
  if value is None:
    return []
  if isinstance(value, (list, tuple)) or value.__class__.__name__ == 'MultiValue':
    return list(value)
  return [value]


def analyze_enhanced_rtstruct(dataset):
  """Return an analysis compatible with dRTImporterPlugin, or None.

  Require the versioned DynamicPET private temporal map.  Static Enhanced PET
  RTSTRUCT objects can also reference multiple frame numbers and must never be
  misidentified as dynamic based only on ReferencedFrameNumber.
  """
  if (str(getattr(dataset, 'SOPClassUID', '')) != RTSTRUCT_SOP_CLASS_UID
      or str(getattr(dataset, 'Modality', '')) != 'RTSTRUCT'):
    return None
  try:
    block = dataset.private_block(0x0011, 'SlicerDynamicPET')
    if block is None:
      return None
    payload = json.loads(str(dataset[block.get_tag(0x01)].value))
  except (KeyError, ValueError, TypeError, AttributeError):
    return None

  if not isinstance(payload, dict):
    return None
  if (type(payload.get('schemaVersion')) is not int
      or payload['schemaVersion'] != 1
      or payload.get('spatialMapping') != 'SINGLE_SOP_ENHANCED_FRAMES'):
    return None
  uid = payload.get('sopInstanceUID')
  sources = payload.get('sourceFrameNumbers')
  if (not isinstance(uid, str) or not uid
      or not isinstance(sources, list) or len(sources) < 2):
    return None

  reverse = {}
  plane_count = None
  for temporal_index, source_frame_numbers in enumerate(sources, start=1):
    if not isinstance(source_frame_numbers, list):
      return None
    if plane_count is None:
      plane_count = len(source_frame_numbers)
    if (plane_count < 2 or len(source_frame_numbers) != plane_count):
      return None
    for k, frame_number in enumerate(source_frame_numbers):
      if (not isinstance(frame_number, Integral) or isinstance(frame_number, bool) or frame_number < 1
          or frame_number in reverse):
        return None
      reverse[frame_number] = (temporal_index, k)

  try:
    reference_frame = dataset.ReferencedFrameOfReferenceSequence[0]
    frame_of_reference_uid = str(reference_frame.FrameOfReferenceUID)
    series = (reference_frame.RTReferencedStudySequence[0]
              .RTReferencedSeriesSequence[0])
    series_uid = str(series.SeriesInstanceUID)
    references = list(series.ContourImageSequence)
    if (len(references) != 1
        or str(references[0].ReferencedSOPInstanceUID) != uid
        or str(references[0].ReferencedSOPClassUID)
        != ENHANCED_PET_SOP_CLASS_UID
        or getattr(references[0], 'ReferencedFrameNumber', None) is not None
        or not frame_of_reference_uid or not series_uid):
      return None
  except (AttributeError, IndexError, TypeError, ValueError):
    return None

  frame_counts = {str(t): 0 for t in range(1, len(sources) + 1)}
  total_contours = 0
  for roi in getattr(dataset, 'ROIContourSequence', []):
    for contour in getattr(roi, 'ContourSequence', []):
      total_contours += 1
      if str(getattr(contour, 'ContourGeometricType', '')) != 'CLOSED_PLANAR':
        return None
      temporal_planes = {}
      contour_refs = list(getattr(contour, 'ContourImageSequence', []))
      if not contour_refs:
        return None
      for reference in contour_refs:
        if (str(getattr(reference, 'ReferencedSOPInstanceUID', '')) != uid
            or str(getattr(reference, 'ReferencedSOPClassUID', ''))
            != ENHANCED_PET_SOP_CLASS_UID):
          return None
        try:
          frame_numbers = _numbers(reference.ReferencedFrameNumber)
        except AttributeError:
          return None
        if not frame_numbers:
          return None
        for frame_number in frame_numbers:
          if not isinstance(frame_number, Integral) or isinstance(frame_number, bool) or frame_number not in reverse:
            return None
          temporal_index, k = reverse[frame_number]
          # A planar contour belongs to one K plane.  The same contour may
          # represent identical geometry at multiple temporal states.
          if temporal_index in temporal_planes:
            return None
          temporal_planes[temporal_index] = k
      if len(set(temporal_planes.values())) != 1:
        return None
      for temporal_index in temporal_planes:
        frame_counts[str(temporal_index)] += 1

  if total_contours == 0:
    return None
  return {
    'convention': 'EnhancedPrivateTemporalMapV1',
    'frameNumbers': list(range(1, len(sources) + 1)),
    'frameToSOPInstanceUIDs': {
      str(t): [uid] for t in range(1, len(sources) + 1)},
    'contourCountByFrame': frame_counts,
    'contourCount': total_contours,
    'referencedSeriesInstanceUID': series_uid,
    'referencedFrameOfReferenceUID': frame_of_reference_uid,
    'enhancedSourceSOPInstanceUID': uid,
    'enhancedSourceFrameNumbers': sources,
    'sourceFrameNumberToTemporalFrame': {
      number: temporal_index for number, (temporal_index, _) in reverse.items()},
  }


def matches_enhanced_pet_sequence(sequence_node, analysis):
  """Match the *full* temporal/spatial map, not repeating SOP UID sets."""
  if sequence_node is None:
    return False
  expected_numbers = analysis.get('enhancedSourceFrameNumbers')
  expected_uid = analysis.get('enhancedSourceSOPInstanceUID')
  if not expected_numbers or not expected_uid:
    return False
  if sequence_node.GetNumberOfDataNodes() != len(expected_numbers):
    return False
  try:
    record = json.loads(sequence_node.GetAttribute('dPET.DICOM.FrameReferences') or '')
    if int(record.get('schemaVersion', 0)) < 2:
      return False
    if str(record.get('seriesInstanceUID') or '') != str(analysis['referencedSeriesInstanceUID']):
      return False
    if str(record.get('frameOfReferenceUID') or '') != str(analysis['referencedFrameOfReferenceUID']):
      return False
    frames = record.get('frames') or []
    if len(frames) != len(expected_numbers):
      return False
    for i, frame in enumerate(frames):
      if frame.get('spatialMapping') != 'SINGLE_SOP_ENHANCED_FRAMES':
        return False
      instances = frame.get('instances') or []
      if (len(instances) != 1
          or str(instances[0].get('sopInstanceUID') or '') != expected_uid
          or str(instances[0].get('sopClassUID') or '') != ENHANCED_PET_SOP_CLASS_UID
          or frame.get('sourceFrameNumbers') != expected_numbers[i]):
        return False
    return True
  except (KeyError, ValueError, TypeError, AttributeError):
    return False
