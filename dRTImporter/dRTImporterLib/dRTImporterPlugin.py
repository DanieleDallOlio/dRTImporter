import hashlib
import json
import logging
import struct

import qt
import slicer
import vtk

from DICOMLib import DICOMLoadable, DICOMPlugin


class dRTImporterPluginClass(DICOMPlugin):
  """DICOM importer for temporal RT Structure Sets used with dynamic PET.

  A supported object stores multiple temporal contour states in one RT
  Structure Set. Referenced Frame Number (0008,1160) is interpreted as the
  one-based temporal PET frame index only when each referenced PET SOP
  Instance UID belongs to exactly one temporal frame. A temporal frame may
  contain one source SOP (3D-per-frame PET) or multiple source SOPs
  (2D-slices-per-frame PET).

  Import is one-pass: pydicom parses the RTSTRUCT once and the plugin builds
  Slicer Planar contours directly. It never writes one temporary RTSTRUCT per
  frame and never invokes the full SlicerRT DICOM loader in a temporal loop.

  Detection is intentionally strict. A static RTSTRUCT must not be claimed by
  this plugin merely because it references several spatial frames of a
  multi-frame image.
  """

  RTSTRUCT_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.481.3"

  def __init__(self):
    super().__init__()
    self.loadType = "dRTImporter"
    self.tags['modality'] = "0008,0060"
    self.tags['sopClassUID'] = "0008,0016"
    self.tags['seriesDescription'] = "0008,103E"
    self.tags['seriesNumber'] = "0020,0011"

  @staticmethod
  def _asIntList(value):
    if value is None:
      return []
    if isinstance(value, (list, tuple)) or value.__class__.__name__ == 'MultiValue':
      return [int(item) for item in value]
    return [int(value)]

  @staticmethod
  def _referencedSeriesInstanceUID(dataset):
    try:
      return str(
        dataset.ReferencedFrameOfReferenceSequence[0]
        .RTReferencedStudySequence[0]
        .RTReferencedSeriesSequence[0]
        .SeriesInstanceUID)
    except Exception:
      return None

  def _analyzeDynamicConvention(self, dataset):
    """Return temporal-reference information, or None if unsupported.

    The accepted convention has a consecutive one-based temporal frame index.
    Every referenced PET SOP Instance UID must belong to exactly one temporal
    frame, while each temporal frame may reference one SOP (3D-per-frame PET)
    or many SOPs (2D-slices-per-frame PET). The complete mapping may be
    declared by the referenced-series ContourImageSequence; contour references
    provide the actual ROI state assignments.
    """
    if (str(getattr(dataset, 'SOPClassUID', '')) != self.RTSTRUCT_SOP_CLASS_UID
        or str(getattr(dataset, 'Modality', '')) != 'RTSTRUCT'):
      return None

    frameToSOPInstanceUIDs = {}
    sopInstanceUIDToFrames = {}
    contourCountByFrame = {}
    contourCount = 0

    def addReference(frameNumber, sopInstanceUID):
      if frameNumber < 1 or not sopInstanceUID:
        return False
      frameToSOPInstanceUIDs.setdefault(frameNumber, set()).add(sopInstanceUID)
      sopInstanceUIDToFrames.setdefault(sopInstanceUID, set()).add(frameNumber)
      return True

    # Optional complete temporal declaration at referenced-series level.
    try:
      referencedSeries = (
        dataset.ReferencedFrameOfReferenceSequence[0]
        .RTReferencedStudySequence[0]
        .RTReferencedSeriesSequence[0])
      for imageReference in getattr(referencedSeries, 'ContourImageSequence', []):
        sopInstanceUID = str(
          getattr(imageReference, 'ReferencedSOPInstanceUID', '') or '')
        for frameNumber in self._asIntList(
            getattr(imageReference, 'ReferencedFrameNumber', None)):
          if not addReference(frameNumber, sopInstanceUID):
            return None
    except Exception:
      pass

    # Contours provide the actual temporal state assignment.
    for roiContour in getattr(dataset, 'ROIContourSequence', []):
      for contour in getattr(roiContour, 'ContourSequence', []):
        contourCount += 1
        contourFrames = set()
        for imageReference in getattr(contour, 'ContourImageSequence', []):
          sopInstanceUID = str(
            getattr(imageReference, 'ReferencedSOPInstanceUID', '') or '')
          frameNumbers = self._asIntList(
            getattr(imageReference, 'ReferencedFrameNumber', None))
          if not sopInstanceUID or not frameNumbers:
            continue
          for frameNumber in frameNumbers:
            if not addReference(frameNumber, sopInstanceUID):
              return None
            contourFrames.add(frameNumber)

        for frameNumber in contourFrames:
          contourCountByFrame[frameNumber] = contourCountByFrame.get(frameNumber, 0) + 1

    frameNumbers = sorted(frameToSOPInstanceUIDs)
    if len(frameNumbers) < 2:
      return None
    if frameNumbers != list(range(1, frameNumbers[-1] + 1)):
      return None
    # The same SOP instance must never be assigned to two temporal states.
    # This is what rejects ordinary enhanced/static RTSTRUCT spatial frame
    # references, while allowing multiple classic 2D PET slices in one state.
    if any(len(frames) != 1 for frames in sopInstanceUIDToFrames.values()):
      return None
    if any(len(uids) < 1 for uids in frameToSOPInstanceUIDs.values()):
      return None

    return {
      'frameNumbers': frameNumbers,
      'frameToSOPInstanceUIDs': {
        str(frame): sorted(frameToSOPInstanceUIDs[frame])
        for frame in frameNumbers
      },
      'contourCountByFrame': {
        str(frame): contourCountByFrame.get(frame, 0)
        for frame in frameNumbers
      },
      'contourCount': contourCount,
      'referencedSeriesInstanceUID': self._referencedSeriesInstanceUID(dataset),
    }

  def examineForImport(self, fileLists):
    return self.examine(fileLists)

  def examine(self, fileLists):
    """Offer only dynamic RTSTRUCT loadables; static RTSTRUCT returns []."""
    import pydicom

    loadables = []
    examinedPaths = set()
    for files in fileLists:
      for filePath in files:
        if filePath in examinedPaths:
          continue
        examinedPaths.add(filePath)
        try:
          dataset = pydicom.dcmread(
            filePath, stop_before_pixels=True, force=True)
          analysis = self._analyzeDynamicConvention(dataset)
          if analysis is None:
            continue

          frameNumbers = analysis['frameNumbers']
          displayName = self._displayBaseName(dataset)

          loadable = DICOMLoadable()
          loadable.files = [filePath]
          loadable.name = (
            f"{displayName} [Dynamic RTSTRUCT] "
            f"({len(frameNumbers)} frames)")
          loadable.tooltip = (
            "Load as a dynamic segmentation sequence using the one-based "
            f"temporal Referenced Frame Number convention (1..{frameNumbers[-1]}).")
          loadable.selected = True
          loadable.confidence = 1.0
          loadables.append(loadable)
        except Exception as error:
          logging.debug(f"[dRTImporter] Examine skipped {filePath}: {error}")
    return loadables

  @staticmethod
  def _sequenceSourceSeriesInstanceUID(sequenceNode):
    seriesInstanceUID = (
      sequenceNode.GetAttribute('dPET.DICOM.SeriesInstanceUID')
      or sequenceNode.GetAttribute('dPET.SourceSeriesInstanceUID')
      or sequenceNode.GetAttribute('DICOM.SeriesInstanceUID'))
    if seriesInstanceUID:
      return seriesInstanceUID

    frameFileList = sequenceNode.GetAttribute('MultiVolume.FrameFileList')
    if frameFileList:
      try:
        import pydicom
        firstFile = frameFileList.split(',')[0]
        dataset = pydicom.dcmread(
          firstFile,
          stop_before_pixels=True,
          force=True,
          specific_tags=['SeriesInstanceUID'])
        return str(dataset.SeriesInstanceUID)
      except Exception:
        pass

    firstDataNode = (
      sequenceNode.GetNthDataNode(0)
      if sequenceNode.GetNumberOfDataNodes() else None)
    instanceUIDs = (
      firstDataNode.GetAttribute('DICOM.instanceUIDs')
      if firstDataNode else None)
    if not instanceUIDs:
      return None
    try:
      import pydicom
      firstFile = slicer.dicomDatabase.fileForInstance(instanceUIDs.split()[0])
      dataset = pydicom.dcmread(
        firstFile,
        stop_before_pixels=True,
        force=True,
        specific_tags=['SeriesInstanceUID'])
      return str(dataset.SeriesInstanceUID)
    except Exception:
      return None

  @staticmethod
  def _sequenceFrameSOPInstanceUIDSets(sequenceNode):
    """Resolve source SOP UID sets for each temporal PET sequence item.

    Priority:
      1. compact provenance persisted by dPETImporter;
      2. DICOM.instanceUIDs already stored on each sequence data node;
      3. legacy dPET.SourceDICOMFiles as a final fallback.
    """
    import os

    provenanceJson = sequenceNode.GetAttribute('dPET.DICOM.FrameReferences')
    if provenanceJson:
      try:
        provenance = json.loads(provenanceJson)
        frames = provenance.get('frames') or []
        if len(frames) == sequenceNode.GetNumberOfDataNodes():
          result = []
          for frame in frames:
            uids = {
              str(item.get('sopInstanceUID') or '')
              for item in (frame.get('instances') or [])
              if str(item.get('sopInstanceUID') or '')
            }
            if not uids:
              raise ValueError('empty provenance frame')
            result.append(uids)
          return result
      except Exception:
        pass

    result = []
    for itemIndex in range(sequenceNode.GetNumberOfDataNodes()):
      dataNode = sequenceNode.GetNthDataNode(itemIndex)
      instanceUIDs = set(
        (dataNode.GetAttribute('DICOM.instanceUIDs') or '').split())
      if instanceUIDs:
        result.append(instanceUIDs)
        continue

      sourceFilesJson = dataNode.GetAttribute('dPET.SourceDICOMFiles')
      if sourceFilesJson:
        try:
          import pydicom
          sourceFiles = json.loads(sourceFilesJson)
          uids = set()
          for filePath in sourceFiles:
            if not os.path.isfile(filePath):
              continue
            dataset = pydicom.dcmread(
              filePath, stop_before_pixels=True, force=True,
              specific_tags=['SOPInstanceUID'])
            uids.add(str(dataset.SOPInstanceUID))
          if uids:
            result.append(uids)
            continue
        except Exception:
          pass

      return None

    return result

  def _findCompatibleVolumeBrowser(self, analysis):
    """Find the uniquely referenced PET/volume browser without guessing."""
    referencedSeriesUID = analysis['referencedSeriesInstanceUID']
    frameNumbers = analysis['frameNumbers']
    expectedSets = [
      frozenset(analysis['frameToSOPInstanceUIDs'][str(frame)])
      for frame in frameNumbers
    ]
    expectedCollection = set(expectedSets)
    if len(expectedCollection) != len(expectedSets):
      return None

    exactMatches = []
    for browserNode in slicer.util.getNodesByClass('vtkMRMLSequenceBrowserNode'):
      masterSequence = browserNode.GetMasterSequenceNode()
      if masterSequence is None:
        continue
      if masterSequence.GetNumberOfDataNodes() != len(frameNumbers):
        continue
      firstDataNode = masterSequence.GetNthDataNode(0)
      if firstDataNode is None or not firstDataNode.IsA('vtkMRMLScalarVolumeNode'):
        continue

      if (referencedSeriesUID
          and self._sequenceSourceSeriesInstanceUID(masterSequence)
          != referencedSeriesUID):
        continue

      sequenceSets = self._sequenceFrameSOPInstanceUIDSets(masterSequence)
      if sequenceSets is None:
        continue
      sequenceCollection = {frozenset(uids) for uids in sequenceSets}
      if (len(sequenceCollection) != len(sequenceSets)
          or sequenceCollection != expectedCollection):
        continue
      exactMatches.append(browserNode)

    if len(exactMatches) == 1:
      return exactMatches[0]
    if len(exactMatches) > 1:
      logging.warning(
        '[dRTImporter] More than one volume browser matches the referenced '
        'PET Series/frame SOP UID sets; creating a standalone browser.')
    else:
      logging.info(
        '[dRTImporter] Referenced PET sequence is not currently available; '
        'creating a standalone dynamic segmentation browser.')
    return None

  @staticmethod
  def _displayBaseName(dataset):
    """Return a stable, human-readable name for the imported dRT object."""
    for attributeName in ('SeriesDescription', 'StructureSetName', 'StructureSetLabel'):
      value = str(getattr(dataset, attributeName, '') or '').strip()
      if value:
        return value
    return 'Dynamic RTSTRUCT'

  @staticmethod
  def _roiDefinitions(dataset):
    """Return ordered DICOM ROI definitions with stable IDs and display colors."""
    colors = {}
    for roiContour in getattr(dataset, 'ROIContourSequence', []):
      try:
        roiNumber = int(roiContour.ReferencedROINumber)
      except Exception:
        continue
      color = getattr(roiContour, 'ROIDisplayColor', None)
      if color is not None and len(color) >= 3:
        colors[roiNumber] = tuple(
          max(0.0, min(1.0, float(component) / 255.0))
          for component in color[:3])

    definitions = []
    for roiItem in getattr(dataset, 'StructureSetROISequence', []):
      roiNumber = int(roiItem.ROINumber)
      definitions.append({
        'number': roiNumber,
        'id': f'dRT_ROI_{roiNumber}',
        'name': str(getattr(roiItem, 'ROIName', '') or f'ROI {roiNumber}'),
        'color': colors.get(roiNumber, (0.5, 0.5, 0.5)),
        'algorithm': str(getattr(roiItem, 'ROIGenerationAlgorithm', '') or ''),
      })
    return definitions

  def _contoursByROIAndFrame(self, dataset, frameNumbers):
    """Index DICOM contours once without duplicating their coordinate arrays."""
    frameSet = set(int(value) for value in frameNumbers)
    indexed = {}
    contourRecords = []

    for roiContour in getattr(dataset, 'ROIContourSequence', []):
      try:
        roiNumber = int(roiContour.ReferencedROINumber)
      except Exception:
        continue

      for contour in getattr(roiContour, 'ContourSequence', []):
        geometricType = str(getattr(contour, 'ContourGeometricType', '') or '')
        if geometricType != 'CLOSED_PLANAR':
          logging.warning(
            '[dRTImporter] Skipping unsupported contour type %s in ROI %s',
            geometricType, roiNumber)
          continue

        values = [float(value) for value in getattr(contour, 'ContourData', [])]
        if len(values) < 9 or len(values) % 3 != 0:
          logging.warning(
            '[dRTImporter] Skipping malformed contour in ROI %s', roiNumber)
          continue
        pointsLPS = [
          (values[index], values[index + 1], values[index + 2])
          for index in range(0, len(values), 3)
        ]
        # Some producers repeat the first vertex explicitly. vtkPolyLine is
        # closed below, so remove that duplicate to avoid a zero-length edge.
        if len(pointsLPS) > 3:
          first = pointsLPS[0]
          last = pointsLPS[-1]
          if sum((first[i] - last[i]) ** 2 for i in range(3)) <= 1e-12:
            pointsLPS = pointsLPS[:-1]

        contourFrames = set()
        for imageReference in getattr(contour, 'ContourImageSequence', []):
          for frameNumber in self._asIntList(
              getattr(imageReference, 'ReferencedFrameNumber', None)):
            if frameNumber in frameSet:
              contourFrames.add(frameNumber)
        if not contourFrames:
          continue

        recordIndex = len(contourRecords)
        contourRecords.append(pointsLPS)
        for frameNumber in contourFrames:
          indexed.setdefault(roiNumber, {}).setdefault(frameNumber, []).append(
            recordIndex)

    return indexed, contourRecords

  @staticmethod
  def _contourStateSignature(contourRecordIndices, contourRecords):
    """Stable geometry digest for one ROI temporal state."""
    digest = hashlib.blake2b(digest_size=16)
    for recordIndex in contourRecordIndices:
      contourPoints = contourRecords[recordIndex]
      digest.update(struct.pack('<I', len(contourPoints)))
      for point in contourPoints:
        digest.update(struct.pack(
          '<3d', float(point[0]), float(point[1]), float(point[2])))
    return digest.hexdigest()

  @staticmethod
  def _planarPolyData(contourRecordIndices, contourRecords):
    """Build SlicerRT-compatible planar-contour vtkPolyData from DICOM LPS."""
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()

    for recordIndex in contourRecordIndices:
      contourPoints = contourRecords[recordIndex]
      if len(contourPoints) < 3:
        continue
      pointIds = []
      for lps in contourPoints:
        # DICOM patient coordinates are LPS; Slicer world coordinates are RAS.
        pointIds.append(points.InsertNextPoint(
          -float(lps[0]), -float(lps[1]), float(lps[2])))

      # SlicerRT's planar contour conversion rules operate on closed line cells.
      # Repeat the first point ID to make closure explicit.
      polyLine = vtk.vtkPolyLine()
      polyLine.GetPointIds().SetNumberOfIds(len(pointIds) + 1)
      for index, pointId in enumerate(pointIds):
        polyLine.GetPointIds().SetId(index, pointId)
      polyLine.GetPointIds().SetId(len(pointIds), pointIds[0])
      lines.InsertNextCell(polyLine)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.SetLines(lines)
    # SlicerRT conversion rules call GetCell()/GetNumberOfCells(); explicitly
    # build the cell links instead of relying on lazy VTK construction.
    polyData.BuildCells()
    polyData.Modified()
    return polyData

  @staticmethod
  def _polyDataSignature(polyData):
    """Exact-enough deterministic geometry fingerprint for diagnostics.

    Includes point coordinates and ordered cell connectivity.  This is used to
    verify that ROI identity remains distinct across DICOM parsing, Slicer
    representation conversion, and sequence storage.
    """
    digest = hashlib.blake2b(digest_size=16)
    if polyData is None:
      digest.update(b'NONE')
      return digest.hexdigest()
    digest.update(struct.pack(
      '<QQ', int(polyData.GetNumberOfPoints()), int(polyData.GetNumberOfCells())))
    for pointIndex in range(polyData.GetNumberOfPoints()):
      point = polyData.GetPoint(pointIndex)
      digest.update(struct.pack(
        '<3d', float(point[0]), float(point[1]), float(point[2])))
    for cellIndex in range(polyData.GetNumberOfCells()):
      cell = polyData.GetCell(cellIndex)
      ids = cell.GetPointIds() if cell is not None else None
      count = ids.GetNumberOfIds() if ids is not None else 0
      digest.update(struct.pack('<I', int(count)))
      for item in range(count):
        digest.update(struct.pack('<q', int(ids.GetId(item))))
    return digest.hexdigest()

  @staticmethod
  def _addClosedSurfaceRepresentations(
      segmentationNode,
      planarName,
      stateSignatureBySegment=None,
      closedSurfaceCache=None):
    """Create Closed surface using Slicer's normal segmentation conversion.

    The previous implementation called CreateRepresentationForOneSegment and
    cached its returned vtkPolyData across temporal states.  Although each
    attached object was deep-copied, that extra cache/one-segment path is not
    needed and makes ROI identity harder to audit.  Here we deliberately mirror
    a conventional SlicerRT segmentation: all Planar contours are present first,
    then vtkSegmentation converts the complete segmentation to Closed surface.

    Planar contours remains the source representation.  The optional cache
    arguments are retained for call compatibility but are intentionally unused.
    """
    import vtkSegmentationCorePython as vtkSegmentationCore

    segmentation = segmentationNode.GetSegmentation()
    closedName = (
      vtkSegmentationCore.vtkSegmentationConverter
      .GetSegmentationClosedSurfaceRepresentationName())

    # Capture planar identity before conversion.
    planarSignatures = {}
    planarBounds = {}
    for segmentId in dRTImporterPluginClass._segmentIDs(segmentation):
      segment = segmentation.GetSegment(segmentId)
      planar = vtk.vtkPolyData.SafeDownCast(segment.GetRepresentation(planarName))
      planarSignatures[segmentId] = dRTImporterPluginClass._polyDataSignature(planar)
      bounds = [0.0] * 6
      if planar is not None and planar.GetNumberOfPoints() > 0:
        planar.GetBounds(bounds)
      planarBounds[segmentId] = [round(float(v), 5) for v in bounds]

    # This is the same high-level conversion path used by normal segmentation
    # workflows.  It lets vtkSegmentation/SlicerRT conversion rules process each
    # segment independently without any cross-ROI Python cache.
    if not segmentationNode.CreateClosedSurfaceRepresentation():
      raise RuntimeError('Failed to create Closed surface representations.')

    closedSignatures = {}
    closedBounds = {}
    for segmentId in dRTImporterPluginClass._segmentIDs(segmentation):
      segment = segmentation.GetSegment(segmentId)
      closed = vtk.vtkPolyData.SafeDownCast(segment.GetRepresentation(closedName))
      if closed is None:
        # Keep representation sets uniform for empty temporal ROI states.
        closed = vtk.vtkPolyData()
        segment.AddRepresentation(closedName, closed)
      closedSignatures[segmentId] = dRTImporterPluginClass._polyDataSignature(closed)
      bounds = [0.0] * 6
      if closed.GetNumberOfPoints() > 0:
        closed.GetBounds(bounds)
      closedBounds[segmentId] = [round(float(v), 5) for v in bounds]

    # A conversion must never make two geometrically distinct planar ROIs become
    # the exact same non-empty closed surface.  Stop immediately if that happens
    # instead of silently displaying the wrong anatomy.
    segmentIds = dRTImporterPluginClass._segmentIDs(segmentation)
    for firstIndex in range(len(segmentIds)):
      firstId = segmentIds[firstIndex]
      for secondIndex in range(firstIndex + 1, len(segmentIds)):
        secondId = segmentIds[secondIndex]
        if planarSignatures[firstId] == planarSignatures[secondId]:
          continue
        firstClosed = segmentation.GetSegment(firstId).GetRepresentation(closedName)
        secondClosed = segmentation.GetSegment(secondId).GetRepresentation(closedName)
        if (firstClosed is not None and secondClosed is not None
            and firstClosed.GetNumberOfPoints() > 0
            and secondClosed.GetNumberOfPoints() > 0
            and closedSignatures[firstId] == closedSignatures[secondId]):
          raise RuntimeError(
            'Closed-surface conversion collapsed distinct ROIs onto identical '
            f'geometry: {firstId} planarBounds={planarBounds[firstId]} and '
            f'{secondId} planarBounds={planarBounds[secondId]}, '
            f'closedBounds={closedBounds[firstId]}.')

    segmentationNode.SetAttribute('dRTImporter.ClosedSurfaceReady', '1')
    segmentationNode.SetAttribute(
      'dRTImporter.PlanarGeometrySignatures',
      json.dumps(planarSignatures, separators=(',', ':')))
    segmentationNode.SetAttribute(
      'dRTImporter.ClosedGeometrySignatures',
      json.dumps(closedSignatures, separators=(',', ':')))
    return closedName

  @staticmethod
  def _segmentIDs(segmentation):
    ids = vtk.vtkStringArray()
    segmentation.GetSegmentIDs(ids)
    return [ids.GetValue(i) for i in range(ids.GetNumberOfValues())]

  @staticmethod
  def _configureClosedSurfaceDisplay(segmentationNode):
    """Prefer Closed surface in both 2D and 3D views, as SlicerRT does."""
    import vtkSegmentationCorePython as vtkSegmentationCore

    if segmentationNode is None:
      return
    segmentationNode.CreateDefaultDisplayNodes()
    displayNode = segmentationNode.GetDisplayNode()
    if displayNode is None:
      return
    closedName = (
      vtkSegmentationCore.vtkSegmentationConverter
      .GetSegmentationClosedSurfaceRepresentationName())
    displayNode.SetPreferredDisplayRepresentationName2D(closedName)
    displayNode.SetPreferredDisplayRepresentationName3D(closedName)
    displayNode.Modified()

  @staticmethod
  def _createFrameSegmentation(
      sourceDataset,
      frameNumber,
      roiDefinitions,
      contoursByROIAndFrame,
      contourRecords,
      polyDataCache,
      closedSurfaceCache,
      referenceVolumeNode=None,
      baseName='Dynamic RTSTRUCT'):
    """Create one temporal segmentation with Planar + Closed representations."""
    import vtkSegmentationCorePython as vtkSegmentationCore

    segmentationNode = slicer.mrmlScene.AddNewNodeByClass(
      'vtkMRMLSegmentationNode', f'{baseName} [frame {frameNumber}]')
    segmentationNode.CreateDefaultDisplayNodes()
    segmentation = segmentationNode.GetSegmentation()
    planarName = (
      vtkSegmentationCore.vtkSegmentationConverter
      .GetSegmentationPlanarContourRepresentationName())
    segmentation.SetSourceRepresentationName(planarName)

    stateSummary = {}
    stateSignatureBySegment = {}
    for definition in roiDefinitions:
      roiNumber = definition['number']
      recordIndices = tuple(
        contoursByROIAndFrame.get(roiNumber, {}).get(frameNumber, []))
      stateSignature = dRTImporterPluginClass._contourStateSignature(
        recordIndices, contourRecords)
      cacheKey = (roiNumber, recordIndices)
      cachedPolyData = polyDataCache.get(cacheKey)
      if cachedPolyData is None:
        cachedPolyData = dRTImporterPluginClass._planarPolyData(
          recordIndices, contourRecords)
        polyDataCache[cacheKey] = cachedPolyData

      # Never attach one vtkPolyData instance to multiple sequence items or
      # segments. Sequence proxies use shallow-copy paths, so sharing mutable
      # representation objects can make derived representations stale or alias.
      polyData = vtk.vtkPolyData()
      polyData.DeepCopy(cachedPolyData)
      polyData.BuildCells()

      segment = vtkSegmentationCore.vtkSegment()
      segment.SetName(definition['name'])
      segment.SetColor(*definition['color'])
      segment.SetTag('dRTImporter.ROINumber', str(roiNumber))
      segment.SetTag('DICOM.ROI.Number', str(roiNumber))
      segment.SetTag('dRTImporter.PlanarStateSignature', stateSignature)
      if definition['algorithm']:
        segment.SetTag('DICOM.ROI.GenerationAlgorithm', definition['algorithm'])
      segment.AddRepresentation(planarName, polyData)
      if not segmentation.AddSegment(segment, definition['id']):
        raise RuntimeError(
          f'Failed to create ROI {roiNumber} in temporal frame {frameNumber}.')

      bounds = [0.0] * 6
      polyData.GetBounds(bounds)
      planarSignature = dRTImporterPluginClass._polyDataSignature(polyData)
      stateSignatureBySegment[definition['id']] = stateSignature
      stateSummary[definition['id']] = {
        'roiNumber': roiNumber,
        'signature': stateSignature,
        'planarSignature': planarSignature,
        'points': int(polyData.GetNumberOfPoints()),
        'cells': int(polyData.GetNumberOfCells()),
        'bounds': [round(float(value), 4) for value in bounds],
      }

    # Direct LPS->RAS planar construction must preserve ROI identity.  If two
    # different DICOM contour states become identical vtkPolyData here, the
    # importer has already corrupted the geometry before any surface conversion.
    summaryIds = list(stateSummary.keys())
    for firstIndex in range(len(summaryIds)):
      firstId = summaryIds[firstIndex]
      for secondIndex in range(firstIndex + 1, len(summaryIds)):
        secondId = summaryIds[secondIndex]
        firstState = stateSummary[firstId]
        secondState = stateSummary[secondId]
        if (firstState['signature'] != secondState['signature']
            and firstState['points'] > 0 and secondState['points'] > 0
            and firstState['planarSignature'] == secondState['planarSignature']):
          raise RuntimeError(
            'Direct planar import collapsed distinct DICOM ROIs onto identical '
            f'geometry in temporal frame {frameNumber}: {firstId} and {secondId}.')

    segmentationNode.SetAttribute('dRTImporter.FrameNumber', str(frameNumber))
    segmentationNode.SetAttribute('dRTImporter.DynamicRTStruct', '1')
    segmentationNode.SetAttribute(
      'dRTImporter.SourceRepresentation', 'Planar contours')
    if referenceVolumeNode is not None:
      segmentationNode.SetReferenceImageGeometryParameterFromVolumeNode(
        referenceVolumeNode)

    # Match normal SlicerRT RTSTRUCT behavior: source remains Planar contours,
    # while Closed surface is derived and already stored in every sequence item.
    dRTImporterPluginClass._addClosedSurfaceRepresentations(
      segmentationNode,
      planarName,
      stateSignatureBySegment,
      closedSurfaceCache)
    segmentationNode.SetAttribute(
      'dRTImporter.GeometrySummary',
      json.dumps(stateSummary, separators=(',', ':')))
    return segmentationNode

  @staticmethod
  def _volumeIndexBySOPInstanceUIDSet(volumeSequence):
    uidSets = dRTImporterPluginClass._sequenceFrameSOPInstanceUIDSets(volumeSequence)
    if uidSets is None:
      return None
    keys = [frozenset(uids) for uids in uidSets]
    if len(set(keys)) != len(keys):
      return None
    return {key: index for index, key in enumerate(keys)}

  def load(self, loadable):
    """Load one dynamic RTSTRUCT in one pass as a planar-contour sequence.

    No temporary DICOM files are created and SlicerRT's full DICOM loader is
    not called per temporal state.  The original contour coordinates are kept
    as Slicer's Planar contours representation.  If the referenced PET is
    already loaded, its frame geometry is attached to each segmentation state
    for later on-demand conversion to binary labelmap.
    """
    import time
    import pydicom

    startTime = time.perf_counter()
    # Ensure SlicerRT's segmentation conversion rules are registered once.
    # This does not invoke the DICOM RT loader or touch the DICOM database.
    if hasattr(slicer.modules, 'dicomrtimportexport'):
      try:
        slicer.modules.dicomrtimportexport.logic()
      except Exception:
        pass
    else:
      logging.warning(
        '[dRTImporter] SlicerRT is not available; Planar contours can be '
        'loaded, but on-demand conversion to binary labelmap may be unavailable.')

    sourcePath = loadable.files[0]
    sourceDataset = pydicom.dcmread(sourcePath, stop_before_pixels=True, force=True)
    analysis = self._analyzeDynamicConvention(sourceDataset)
    if analysis is None:
      logging.info(
        '[dRTImporter] The selected RT Structure Set is not dynamic; ignored.')
      return False

    frameNumbers = analysis['frameNumbers']
    roiDefinitions = self._roiDefinitions(sourceDataset)
    if not roiDefinitions:
      logging.error('[dRTImporter] Dynamic RTSTRUCT contains no ROI definitions.')
      return False

    volumeBrowser = self._findCompatibleVolumeBrowser(analysis)
    volumeSequence = (
      volumeBrowser.GetMasterSequenceNode() if volumeBrowser else None)
    volumeIndexByUIDSet = (
      self._volumeIndexBySOPInstanceUIDSet(volumeSequence)
      if volumeSequence else None)

    baseName = self._displayBaseName(sourceDataset)
    segmentationSequence = slicer.mrmlScene.AddNewNodeByClass(
      'vtkMRMLSequenceNode',
      slicer.mrmlScene.GenerateUniqueName(baseName + ' [dRT]'))
    segmentationSequence.SetIndexName(
      volumeSequence.GetIndexName() if volumeSequence else 'frame')
    segmentationSequence.SetIndexUnit(
      volumeSequence.GetIndexUnit() if volumeSequence else '')
    segmentationSequence.SetIndexType(
      volumeSequence.GetIndexType() if volumeSequence else 0)
    segmentationSequence.SetAttribute('dRTImporter.DynamicRTStruct', '1')
    segmentationSequence.SetAttribute('dRTImporter.DisplayName', baseName)
    segmentationSequence.SetAttribute(
      'dRTImporter.Convention', 'TemporalReferencedFrameNumber')
    segmentationSequence.SetAttribute('dRTImporter.SourceFile', sourcePath)
    segmentationSequence.SetAttribute(
      'dRTImporter.SourceRepresentation', 'Planar contours')
    segmentationSequence.SetAttribute(
      'dRTImporter.FrameToSOPInstanceUIDs',
      json.dumps(analysis['frameToSOPInstanceUIDs'], separators=(',', ':')))
    if analysis['referencedSeriesInstanceUID']:
      segmentationSequence.SetAttribute(
        'dRTImporter.ReferencedSeriesInstanceUID',
        analysis['referencedSeriesInstanceUID'])

    progress = slicer.util.createProgressDialog(
      labelText=f'Reading {baseName}',
      value=0,
      maximum=100,
      windowModality=qt.Qt.WindowModal)

    try:
      progress.value = 5
      progress.labelText = 'Indexing dynamic RTSTRUCT contours'
      slicer.app.processEvents()
      contoursByROIAndFrame, contourRecords = self._contoursByROIAndFrame(
        sourceDataset, frameNumbers)

      # Summarize temporal geometry states before any VTK conversion. This
      # proves whether an edited ROI state is present in the DICOM object.
      roiStateSignatures = {definition['number']: set() for definition in roiDefinitions}
      for definition in roiDefinitions:
        roiNumber = definition['number']
        for frameNumber in frameNumbers:
          recordIndices = tuple(
            contoursByROIAndFrame.get(roiNumber, {}).get(frameNumber, []))
          roiStateSignatures[roiNumber].add(
            self._contourStateSignature(recordIndices, contourRecords))
      logging.info(
        '[dRTImporter] Temporal ROI states from DICOM: %s',
        ', '.join(
          f"ROI {definition['number']} '{definition['name']}'="
          f"{len(roiStateSignatures[definition['number']])} unique state(s)"
          for definition in roiDefinitions))

      if progress.wasCanceled:
        raise RuntimeError('Dynamic RTSTRUCT import was cancelled.')

      polyDataCache = {}
      closedSurfaceCache = None
      for itemIndex, frameNumber in enumerate(frameNumbers):
        progress.value = 10 + int(round(80.0 * (itemIndex + 1) / max(1, len(frameNumbers))))
        progress.labelText = (
          f'Building dRT frame {itemIndex + 1}/{len(frameNumbers)} ' \
          '(including Closed surface)')
        slicer.app.processEvents()
        if progress.wasCanceled:
          raise RuntimeError('Dynamic RTSTRUCT import was cancelled.')

        referenceVolumeNode = None
        if volumeSequence:
          referencedUIDSet = frozenset(
            analysis['frameToSOPInstanceUIDs'][str(frameNumber)])
          if (volumeIndexByUIDSet is None
              or referencedUIDSet not in volumeIndexByUIDSet):
            raise RuntimeError(
              f'PET SOP Instance UID set for dynamic RT frame {frameNumber} '
              'cannot be mapped to the selected PET sequence.')
          volumeItemNumber = volumeIndexByUIDSet[referencedUIDSet]
          indexValue = volumeSequence.GetNthIndexValue(volumeItemNumber)
          referenceVolumeNode = volumeSequence.GetNthDataNode(volumeItemNumber)
        else:
          indexValue = str(frameNumber - 1)

        segmentationNode = self._createFrameSegmentation(
          sourceDataset,
          frameNumber,
          roiDefinitions,
          contoursByROIAndFrame,
          contourRecords,
          polyDataCache,
          closedSurfaceCache,
          referenceVolumeNode,
          baseName)
        segmentationSequence.SetDataNodeAtValue(segmentationNode, indexValue)
        slicer.mrmlScene.RemoveNode(segmentationNode)

      progress.value = 92
      progress.labelText = 'Connecting dynamic segmentation sequence'
      slicer.app.processEvents()

      if volumeBrowser:
        browser = volumeBrowser
        browser.AddSynchronizedSequenceNode(segmentationSequence)
      else:
        browser = slicer.mrmlScene.AddNewNodeByClass(
          'vtkMRMLSequenceBrowserNode',
          slicer.mrmlScene.GenerateUniqueName(baseName + ' [dRT Browser]'))
        browser.SetAndObserveMasterSequenceNodeID(segmentationSequence.GetID())

      browser.SetSaveChanges(segmentationSequence, True)
      # Keep one stable human-readable proxy name while the browser changes frames.
      # Overwriting it from sequence item names on every slider move caused
      # repeated 'sequence' suffixes to accumulate.
      browser.SetOverwriteProxyName(segmentationSequence, False)
      browser.SetAttribute(
        'dRTImporter.SequenceNodeID', segmentationSequence.GetID())

      proxySegmentation = browser.GetProxyNode(segmentationSequence)
      if proxySegmentation:
        proxySegmentation.SetName(baseName)
        proxySegmentation.SetAttribute('dRTImporter.DynamicRTStruct', '1')
        proxySegmentation.SetAttribute(
          'dRTImporter.SequenceNodeID', segmentationSequence.GetID())
        proxySegmentation.SetAttribute(
          'dRTImporter.SourceRepresentation', 'Planar contours')
        proxySegmentation.SetAttribute(
          'dRTImporter.FrameToSOPInstanceUIDs',
          json.dumps(analysis['frameToSOPInstanceUIDs'], separators=(',', ':')))
        if analysis['referencedSeriesInstanceUID']:
          proxySegmentation.SetAttribute(
            'dRTImporter.ReferencedSeriesInstanceUID',
            analysis['referencedSeriesInstanceUID'])
        if volumeBrowser and volumeSequence:
          volumeProxy = volumeBrowser.GetProxyNode(volumeSequence)
          if volumeProxy:
            proxySegmentation.SetReferenceImageGeometryParameterFromVolumeNode(
              volumeProxy)
        self._configureClosedSurfaceDisplay(proxySegmentation)
        try:
          self.addSeriesInSubjectHierarchy(loadable, proxySegmentation)
        except Exception as error:
          logging.warning(
            '[dRTImporter] Subject hierarchy DICOM placement skipped: %s', error)

      if (hasattr(slicer.modules, 'sequences')
          and slicer.modules.sequences.widgetRepresentation()):
        slicer.modules.sequences.widgetRepresentation().setActiveBrowserNode(
          browser)

      progress.value = 100
      progress.labelText = 'Dynamic RTSTRUCT import complete'
      slicer.app.processEvents()
      logging.info(
        '[dRTImporter] Direct Planar+Closed import complete in %.3f s; '
        'frames=%d ROIs=%d contour records=%d unique planar states=%d',
        time.perf_counter() - startTime,
        len(frameNumbers), len(roiDefinitions), len(contourRecords),
        len(polyDataCache))
      logging.info(
        '[dRTImporter] Closed-surface conversion uses standard per-frame segmentation conversion (no cross-ROI cache).')
      return segmentationSequence

    except Exception as error:
      logging.error(f'[dRTImporter] Import failed: {error}')
      import traceback
      traceback.print_exc()
      slicer.mrmlScene.RemoveNode(segmentationSequence)
      return False
    finally:
      progress.close()


def adopt_dynamic_rtstruct_to_pet(
    segmentation_proxy_node_id,
    pet_sequence_node_id,
    pet_browser_node_id):
  """Adopt a standalone dRTImporter sequence into the selected PET browser.

  Returns a Python dictionary suitable for conversion to QVariantMap.
  The operation is conservative: both the referenced PET series and the full
  per-frame SOP Instance UID set must match. A mismatch is reported and the
  segmentation is left untouched.
  """
  try:
    proxyNode = slicer.mrmlScene.GetNodeByID(str(segmentation_proxy_node_id))
    petSequence = slicer.mrmlScene.GetNodeByID(str(pet_sequence_node_id))
    petBrowser = slicer.mrmlScene.GetNodeByID(str(pet_browser_node_id))
    if proxyNode is None or not proxyNode.IsA('vtkMRMLSegmentationNode'):
      raise RuntimeError('The selected segmentation proxy is unavailable.')
    if petSequence is None or not petSequence.IsA('vtkMRMLSequenceNode'):
      raise RuntimeError('The selected PET sequence is unavailable.')
    if petBrowser is None or not petBrowser.IsA('vtkMRMLSequenceBrowserNode'):
      raise RuntimeError('The selected PET sequence browser is unavailable.')

    sequenceNodeID = proxyNode.GetAttribute('dRTImporter.SequenceNodeID')
    if not sequenceNodeID:
      raise RuntimeError('The selected segmentation has no dRTImporter sequence reference.')
    sourceSequence = slicer.mrmlScene.GetNodeByID(sequenceNodeID)
    if sourceSequence is None or not sourceSequence.IsA('vtkMRMLSequenceNode'):
      raise RuntimeError('The imported dynamic segmentation sequence is unavailable.')

    referencedSeriesUID = sourceSequence.GetAttribute(
      'dRTImporter.ReferencedSeriesInstanceUID')
    petSeriesUID = dRTImporterPluginClass._sequenceSourceSeriesInstanceUID(petSequence)
    if (not referencedSeriesUID or not petSeriesUID
        or referencedSeriesUID != petSeriesUID):
      raise RuntimeError(
        'The dynamic RTSTRUCT references a different PET Series Instance UID '
        'than the PET selected in SlicerDynamicPET.')

    frameMapJson = sourceSequence.GetAttribute(
      'dRTImporter.FrameToSOPInstanceUIDs')
    if not frameMapJson:
      # Backward compatibility with the earlier one-SOP-per-frame importer.
      legacyJson = sourceSequence.GetAttribute('dRTImporter.FrameToSOPInstanceUID')
      if legacyJson:
        legacy = json.loads(legacyJson)
        frameMap = {key: [value] for key, value in legacy.items()}
      else:
        raise RuntimeError('The dynamic RTSTRUCT has no frame-to-SOP reference map.')
    else:
      frameMap = json.loads(frameMapJson)

    frameNumbers = sorted(int(value) for value in frameMap.keys())
    referencedSets = [
      frozenset(frameMap[str(frameNumber)]) for frameNumber in frameNumbers
    ]

    petUIDSets = dRTImporterPluginClass._sequenceFrameSOPInstanceUIDSets(petSequence)
    if petUIDSets is None:
      raise RuntimeError(
        'The selected PET sequence does not expose source SOP Instance UID '
        'sets for each temporal frame.')
    petKeys = [frozenset(uids) for uids in petUIDSets]
    if (len(petKeys) != len(referencedSets)
        or set(petKeys) != set(referencedSets)):
      raise RuntimeError(
        'The dynamic RTSTRUCT temporal SOP reference sets do not match the '
        'selected PET sequence.')
    petIndexByUIDSet = {key: index for index, key in enumerate(petKeys)}

    sourceNodeByFrame = {}
    for itemIndex in range(sourceSequence.GetNumberOfDataNodes()):
      dataNode = sourceSequence.GetNthDataNode(itemIndex)
      frameNumberText = dataNode.GetAttribute('dRTImporter.FrameNumber')
      if frameNumberText:
        sourceNodeByFrame[int(frameNumberText)] = dataNode
      elif itemIndex < len(frameNumbers):
        sourceNodeByFrame[frameNumbers[itemIndex]] = dataNode

    if any(frameNumber not in sourceNodeByFrame for frameNumber in frameNumbers):
      raise RuntimeError('The imported segmentation sequence is missing temporal frames.')

    displayName = (
      sourceSequence.GetAttribute('dRTImporter.DisplayName')
      or proxyNode.GetName()
      or 'Dynamic RTSTRUCT')
    alignedSequence = slicer.mrmlScene.AddNewNodeByClass(
      'vtkMRMLSequenceNode',
      slicer.mrmlScene.GenerateUniqueName(displayName + ' [dRT]'))
    alignedSequence.CopySequenceIndex(petSequence)
    for attributeName in (
        'dRTImporter.DynamicRTStruct',
        'dRTImporter.DisplayName',
        'dRTImporter.Convention',
        'dRTImporter.SourceFile',
        'dRTImporter.FrameToSOPInstanceUIDs',
        'dRTImporter.FrameToSOPInstanceUID',
        'dRTImporter.ReferencedSeriesInstanceUID'):
      value = sourceSequence.GetAttribute(attributeName)
      if value is not None:
        alignedSequence.SetAttribute(attributeName, value)

    for frameNumber in frameNumbers:
      referencedUIDSet = frozenset(frameMap[str(frameNumber)])
      petItemIndex = petIndexByUIDSet[referencedUIDSet]
      indexValue = petSequence.GetNthIndexValue(petItemIndex)
      sourceNode = sourceNodeByFrame[frameNumber]
      petFrameNode = petSequence.GetNthDataNode(petItemIndex)
      if petFrameNode is not None:
        sourceNode.SetReferenceImageGeometryParameterFromVolumeNode(petFrameNode)
      # Backward compatibility for sequences imported by an earlier build that
      # stored only Planar contours.
      if sourceNode.GetAttribute('dRTImporter.ClosedSurfaceReady') != '1':
        import vtkSegmentationCorePython as vtkSegmentationCore
        planarName = (
          vtkSegmentationCore.vtkSegmentationConverter
          .GetSegmentationPlanarContourRepresentationName())
        dRTImporterPluginClass._addClosedSurfaceRepresentations(
          sourceNode, planarName)
      alignedSequence.SetDataNodeAtValue(sourceNode, indexValue)

    # Detach the selected proxy from its standalone importer browser while
    # keeping the proxy node itself alive, then reuse it in the PET browser.
    oldBrowsers = []
    for browserNode in slicer.util.getNodesByClass('vtkMRMLSequenceBrowserNode'):
      if browserNode.GetID() == petBrowser.GetID():
        continue
      browserSequence = browserNode.GetSequenceNode(proxyNode)
      if (browserSequence is not None
          and browserSequence.GetID() == sourceSequence.GetID()):
        oldBrowsers.append(browserNode)

    for browserNode in oldBrowsers:
      browserNode.SetSaveChanges(sourceSequence, False)
      proxyNode.RemoveAttribute('proxyNodeCopy')
      browserNode.RemoveSynchronizedSequenceNode(sourceSequence.GetID())
      if (browserNode.GetNumberOfSynchronizedSequenceNodes(True) == 0
          and browserNode.GetAttribute('dRTImporter.SequenceNodeID')):
        slicer.mrmlScene.RemoveNode(browserNode)

    petBrowser.AddProxyNode(proxyNode, alignedSequence, False)
    petBrowser.SetSaveChanges(alignedSequence, True)
    petBrowser.SetOverwriteProxyName(alignedSequence, False)
    proxyNode.SetName(displayName)
    petProxy = petBrowser.GetProxyNode(petSequence)
    if petProxy is not None:
      proxyNode.SetReferenceImageGeometryParameterFromVolumeNode(petProxy)
    proxyNode.SetAttribute('dRTImporter.SequenceNodeID', alignedSequence.GetID())
    dRTImporterPluginClass._configureClosedSurfaceDisplay(proxyNode)
    alignedSequence.SetAttribute('dRTImporter.AdoptedToPET', '1')

    if sourceSequence is not alignedSequence:
      slicer.mrmlScene.RemoveNode(sourceSequence)

    return {
      'ok': True,
      'sequence_node_id': alignedSequence.GetID(),
      'proxy_node_id': proxyNode.GetID(),
    }
  except Exception as error:
    logging.error(f'[dRTImporter] PET adoption failed: {error}')
    return {'ok': False, 'error': str(error)}


__all__ = ['dRTImporterPluginClass', 'adopt_dynamic_rtstruct_to_pet']
