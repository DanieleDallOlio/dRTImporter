import json
import logging

import qt
import slicer
import vtk

from DICOMLib import DICOMLoadable, DICOMPlugin


class dRTImporterPluginClass(DICOMPlugin):
  """DICOM importer for temporal RT Structure Sets used with dynamic PET.

  A supported object stores multiple temporal contour states in one RT
  Structure Set. Referenced Frame Number (0008,1160) is interpreted as the
  one-based temporal PET frame index only when the RTSTRUCT also provides a
  one-to-one mapping between those frame numbers and distinct referenced PET
  SOP Instance UIDs.

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

    The accepted convention has a consecutive one-based temporal frame index
    and a one-to-one mapping between each temporal frame and one distinct PET
    SOP Instance UID. The mapping may be declared by contour references, by a
    complete series-level ContourImageSequence, or by both. Series-level
    references make intentionally empty temporal segmentation states possible.
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
    if any(len(uids) != 1 for uids in frameToSOPInstanceUIDs.values()):
      return None
    if any(len(frames) != 1 for frames in sopInstanceUIDToFrames.values()):
      return None
    if len(frameToSOPInstanceUIDs) != len(sopInstanceUIDToFrames):
      return None

    return {
      'frameNumbers': frameNumbers,
      'frameToSOPInstanceUID': {
        str(frame): next(iter(frameToSOPInstanceUIDs[frame]))
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
          seriesDescription = str(
            getattr(dataset, 'SeriesDescription', '') or 'Dynamic RTSTRUCT')

          loadable = DICOMLoadable()
          loadable.files = [filePath]
          loadable.name = (
            f"{seriesDescription} - dynamic RTSTRUCT "
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
    seriesInstanceUID = sequenceNode.GetAttribute(
      'dPET.SourceSeriesInstanceUID')
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
  def _sequenceFrameSOPInstanceUIDs(sequenceNode):
    """Resolve one source SOP Instance UID per sequence item, in item order."""
    import os
    import pydicom

    result = []
    for itemIndex in range(sequenceNode.GetNumberOfDataNodes()):
      dataNode = sequenceNode.GetNthDataNode(itemIndex)
      uid = None

      sourceFilesJson = dataNode.GetAttribute('dPET.SourceDICOMFiles')
      if sourceFilesJson:
        try:
          sourceFiles = json.loads(sourceFilesJson)
          if len(sourceFiles) == 1 and os.path.isfile(sourceFiles[0]):
            dataset = pydicom.dcmread(
              sourceFiles[0], stop_before_pixels=True, force=True,
              specific_tags=['SOPInstanceUID'])
            uid = str(dataset.SOPInstanceUID)
        except Exception:
          uid = None

      if uid is None:
        instanceUIDs = (dataNode.GetAttribute('DICOM.instanceUIDs') or '').split()
        if len(instanceUIDs) == 1:
          uid = instanceUIDs[0]

      if uid is None:
        return None
      result.append(uid)

    return result

  def _findCompatibleVolumeBrowser(self, analysis):
    """Find the uniquely referenced PET/volume browser without guessing."""
    referencedSeriesUID = analysis['referencedSeriesInstanceUID']
    frameNumbers = analysis['frameNumbers']
    expectedUIDs = [
      analysis['frameToSOPInstanceUID'][str(frame)] for frame in frameNumbers
    ]
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

      sequenceUIDs = self._sequenceFrameSOPInstanceUIDs(masterSequence)
      if sequenceUIDs is None:
        continue
      if (len(sequenceUIDs) != len(expectedUIDs)
          or len(set(sequenceUIDs)) != len(sequenceUIDs)
          or set(sequenceUIDs) != set(expectedUIDs)):
        continue
      exactMatches.append(browserNode)

    if len(exactMatches) == 1:
      return exactMatches[0]
    if len(exactMatches) > 1:
      logging.warning(
        '[dRTImporter] More than one volume browser matches the referenced '
        'PET Series/SOP Instance UID map; creating a standalone browser.')
    else:
      logging.info(
        '[dRTImporter] Referenced PET sequence is not currently available; '
        'creating a standalone dynamic segmentation browser.')
    return None

  @staticmethod
  def _normalizeSegmentIDs(segmentationNode, frameDataset):
    """Replace transient SlicerRT segment IDs with stable DICOM ROI IDs."""
    import vtkSegmentationCorePython as vtkSegmentationCore

    segmentation = segmentationNode.GetSegmentation()
    roiItems = list(getattr(frameDataset, 'StructureSetROISequence', []))
    if not roiItems:
      return

    importedIDs = vtk.vtkStringArray()
    segmentation.GetSegmentIDs(importedIDs)
    idsByName = {}
    for index in range(importedIDs.GetNumberOfValues()):
      segmentID = importedIDs.GetValue(index)
      segment = segmentation.GetSegment(segmentID)
      idsByName.setdefault(segment.GetName(), []).append(segmentID)

    replacements = []
    usedImportedIDs = set()
    for roiItem in roiItems:
      roiNumber = int(roiItem.ROINumber)
      roiName = str(roiItem.ROIName)
      candidates = [
        segmentID for segmentID in idsByName.get(roiName, [])
        if segmentID not in usedImportedIDs
      ]
      if len(candidates) != 1:
        raise RuntimeError(
          f'Cannot uniquely match imported ROI {roiNumber} ({roiName!r}) '
          'to a Slicer segment.')
      sourceID = candidates[0]
      usedImportedIDs.add(sourceID)
      sourceSegment = segmentation.GetSegment(sourceID)
      segmentCopy = vtkSegmentationCore.vtkSegment()
      segmentCopy.DeepCopy(sourceSegment)
      segmentCopy.SetTag('dRTImporter.ROINumber', str(roiNumber))
      replacements.append((f'dRT_ROI_{roiNumber}', segmentCopy))

    # Preserve any unexpected segment rather than silently dropping it.
    for index in range(importedIDs.GetNumberOfValues()):
      sourceID = importedIDs.GetValue(index)
      if sourceID in usedImportedIDs:
        continue
      segmentCopy = vtkSegmentationCore.vtkSegment()
      segmentCopy.DeepCopy(segmentation.GetSegment(sourceID))
      replacements.append((sourceID, segmentCopy))

    segmentation.RemoveAllSegments()
    for stableID, segmentCopy in replacements:
      segmentation.AddSegment(segmentCopy, stableID)

  @staticmethod
  def _ensureAllROIs(segmentationNode, sourceDataset):
    """Ensure every DICOM ROI exists in every temporal segmentation state."""
    import vtkSegmentationCorePython as vtkSegmentationCore

    segmentation = segmentationNode.GetSegmentation()
    colorsByROINumber = {}
    for roiContour in getattr(sourceDataset, 'ROIContourSequence', []):
      try:
        roiNumber = int(roiContour.ReferencedROINumber)
      except Exception:
        continue
      color = getattr(roiContour, 'ROIDisplayColor', None)
      if color is not None and len(color) >= 3:
        colorsByROINumber[roiNumber] = tuple(
          max(0.0, min(1.0, float(component) / 255.0))
          for component in color[:3])

    for roiItem in getattr(sourceDataset, 'StructureSetROISequence', []):
      roiNumber = int(roiItem.ROINumber)
      stableID = f'dRT_ROI_{roiNumber}'
      if segmentation.GetSegment(stableID) is not None:
        continue
      segment = vtkSegmentationCore.vtkSegment()
      segment.SetName(str(roiItem.ROIName))
      color = colorsByROINumber.get(roiNumber, (0.5, 0.5, 0.5))
      segment.SetColor(*color)
      segment.SetTag('dRTImporter.ROINumber', str(roiNumber))
      segmentation.AddSegment(segment, stableID)

  @staticmethod
  def _createEmptyFrameSegmentation(sourceDataset):
    segmentationNode = slicer.mrmlScene.AddNewNodeByClass(
      'vtkMRMLSegmentationNode', 'dRTImporter empty temporal frame')
    segmentationNode.CreateDefaultDisplayNodes()
    dRTImporterPluginClass._ensureAllROIs(segmentationNode, sourceDataset)
    return segmentationNode

  @staticmethod
  def _volumeIndexBySOPInstanceUID(volumeSequence):
    uids = dRTImporterPluginClass._sequenceFrameSOPInstanceUIDs(volumeSequence)
    if uids is None:
      return None
    if len(set(uids)) != len(uids):
      return None
    return {uid: index for index, uid in enumerate(uids)}

  @staticmethod
  def _sanitizeDateTimeValues(dataset):
    """Clear invalid anonymized DA/TM values before strict SlicerRT parsing."""
    from datetime import datetime

    for element in dataset.iterall():
      value = str(element.value or '')
      if not value or element.VR not in ('DA', 'TM'):
        continue
      try:
        if element.VR == 'DA':
          datetime.strptime(value, '%Y%m%d')
        else:
          mainValue = value.split('.')[0]
          datetime.strptime(mainValue.ljust(6, '0'), '%H%M%S')
      except Exception:
        element.value = ''

  def _datasetForFrame(self, sourceDataset, frameNumber):
    """Create one conventional, static RTSTRUCT for SlicerRT conversion."""
    import copy

    from pydicom.sequence import Sequence
    from pydicom.uid import generate_uid

    dataset = copy.deepcopy(sourceDataset)
    retainedROIContours = []
    retainedROINumbers = set()

    for roiContour in getattr(dataset, 'ROIContourSequence', []):
      retainedContours = []
      for contour in getattr(roiContour, 'ContourSequence', []):
        retainContour = False
        for imageReference in getattr(
            contour, 'ContourImageSequence', []):
          if frameNumber in self._asIntList(
              getattr(imageReference, 'ReferencedFrameNumber', None)):
            retainContour = True
            # The temporary object contains only one temporal state, so the
            # temporal marker is no longer needed by the standard reader.
            del imageReference.ReferencedFrameNumber
        if retainContour:
          # Keep only image references that belong to this temporal state.
          retainedReferences = []
          for imageReference in getattr(contour, 'ContourImageSequence', []):
            frameValues = self._asIntList(
              getattr(imageReference, 'ReferencedFrameNumber', None))
            if not frameValues:
              # Matching references have already had the temporal marker
              # removed above.
              retainedReferences.append(imageReference)
            elif frameNumber in frameValues:
              del imageReference.ReferencedFrameNumber
              retainedReferences.append(imageReference)
          contour.ContourImageSequence = Sequence(retainedReferences)
          retainedContours.append(contour)

      if retainedContours:
        roiContour.ContourSequence = Sequence(retainedContours)
        retainedROIContours.append(roiContour)
        retainedROINumbers.add(int(roiContour.ReferencedROINumber))

    dataset.ROIContourSequence = Sequence(retainedROIContours)
    dataset.StructureSetROISequence = Sequence([
      item for item in getattr(dataset, 'StructureSetROISequence', [])
      if int(item.ROINumber) in retainedROINumbers
    ])
    if hasattr(dataset, 'RTROIObservationsSequence'):
      dataset.RTROIObservationsSequence = Sequence([
        item for item in dataset.RTROIObservationsSequence
        if (hasattr(item, 'ReferencedROINumber')
            and int(item.ReferencedROINumber) in retainedROINumbers)
      ])

    # If the dynamic object carries a complete temporal map at referenced-
    # series level, reduce it to this one state for standard SlicerRT import.
    try:
      referencedSeries = (
        dataset.ReferencedFrameOfReferenceSequence[0]
        .RTReferencedStudySequence[0]
        .RTReferencedSeriesSequence[0])
      retainedSeriesReferences = []
      for imageReference in getattr(referencedSeries, 'ContourImageSequence', []):
        frameValues = self._asIntList(
          getattr(imageReference, 'ReferencedFrameNumber', None))
        if frameNumber not in frameValues:
          continue
        del imageReference.ReferencedFrameNumber
        retainedSeriesReferences.append(imageReference)
      if retainedSeriesReferences:
        referencedSeries.ContourImageSequence = Sequence(retainedSeriesReferences)
    except Exception:
      pass

    temporarySOPInstanceUID = generate_uid()
    dataset.SOPInstanceUID = temporarySOPInstanceUID
    dataset.SeriesInstanceUID = generate_uid()
    if getattr(dataset, 'file_meta', None) is not None:
      dataset.file_meta.MediaStorageSOPInstanceUID = temporarySOPInstanceUID
    dataset.SeriesDescription = (
      f"{getattr(sourceDataset, 'SeriesDescription', 'Dynamic RTSTRUCT')} "
      f"frame {frameNumber}")
    self._sanitizeDateTimeValues(dataset)
    return dataset

  @staticmethod
  def _loadTemporaryRTStruct(filePath):
    """Load one filtered RTSTRUCT through the SlicerRT reader."""
    if not hasattr(slicer.modules, 'dicomrtimportexport'):
      raise RuntimeError(
        'SlicerRT is required. Install the SlicerRT extension and restart '
        '3D Slicer.')

    previousIDs = {
      node.GetID()
      for node in slicer.util.getNodesByClass('vtkMRMLSegmentationNode')
    }

    qtLoadable = slicer.qSlicerDICOMLoadable()
    qtLoadable.files = [filePath]
    qtLoadable.name = 'dRTImporter temporary frame'
    qtLoadable.tooltip = qtLoadable.name
    qtLoadable.selected = True
    qtLoadable.confidence = 1.0

    vtkLoadable = slicer.vtkSlicerDICOMLoadable()
    qtLoadable.copyToVtkLoadable(vtkLoadable)
    success = slicer.modules.dicomrtimportexport.logic().LoadDicomRT(
      vtkLoadable)
    if not success:
      raise RuntimeError(
        f'SlicerRT failed to load temporary RTSTRUCT: {filePath}')

    newNodes = [
      node
      for node in slicer.util.getNodesByClass('vtkMRMLSegmentationNode')
      if node.GetID() not in previousIDs
    ]
    if not newNodes:
      raise RuntimeError(
        'SlicerRT reported success but created no segmentation node.')

    selectedNode = max(
      newNodes,
      key=lambda node: node.GetSegmentation().GetNumberOfSegments())
    for node in newNodes:
      if node is not selectedNode:
        slicer.mrmlScene.RemoveNode(node)
    return selectedNode

  def load(self, loadable):
    """Load the dynamic RTSTRUCT as a segmentation sequence."""
    import os
    import tempfile

    import pydicom

    sourcePath = loadable.files[0]
    sourceDataset = pydicom.dcmread(sourcePath, force=True)
    analysis = self._analyzeDynamicConvention(sourceDataset)
    if analysis is None:
      # Defensive second check: a static RTSTRUCT must never be imported.
      logging.info(
        '[dRTImporter] The selected RT Structure Set is not dynamic; ignored.')
      return False

    frameNumbers = analysis['frameNumbers']
    volumeBrowser = self._findCompatibleVolumeBrowser(analysis)
    volumeSequence = (
      volumeBrowser.GetMasterSequenceNode() if volumeBrowser else None)

    baseName = str(
      getattr(sourceDataset, 'SeriesDescription', '')
      or 'Dynamic RTSTRUCT')
    segmentationSequence = slicer.mrmlScene.AddNewNodeByClass(
      'vtkMRMLSequenceNode',
      slicer.mrmlScene.GenerateUniqueName(
        baseName + ' segmentation sequence'))
    segmentationSequence.SetIndexName(
      volumeSequence.GetIndexName() if volumeSequence else 'frame')
    segmentationSequence.SetIndexUnit(
      volumeSequence.GetIndexUnit() if volumeSequence else '')
    segmentationSequence.SetIndexType(
      volumeSequence.GetIndexType() if volumeSequence else 0)
    segmentationSequence.SetAttribute('dRTImporter.DynamicRTStruct', '1')
    segmentationSequence.SetAttribute(
      'dRTImporter.Convention', 'TemporalReferencedFrameNumber')
    segmentationSequence.SetAttribute('dRTImporter.SourceFile', sourcePath)
    segmentationSequence.SetAttribute(
      'dRTImporter.FrameToSOPInstanceUID',
      json.dumps(analysis['frameToSOPInstanceUID']))
    if analysis['referencedSeriesInstanceUID']:
      segmentationSequence.SetAttribute(
        'dRTImporter.ReferencedSeriesInstanceUID',
        analysis['referencedSeriesInstanceUID'])

    progress = slicer.util.createProgressDialog(
      labelText=f'Loading {baseName}',
      value=0,
      maximum=len(frameNumbers),
      windowModality=qt.Qt.WindowModal)

    try:
      with tempfile.TemporaryDirectory(
          prefix='dRTImporter_') as temporaryDirectory:
        for itemIndex, frameNumber in enumerate(frameNumbers):
          progress.value = itemIndex
          progress.labelText = (
            f'Loading RTSTRUCT frame {frameNumber} '
            f'({itemIndex + 1}/{len(frameNumbers)})')
          slicer.app.processEvents()
          if progress.wasCanceled:
            raise RuntimeError('Dynamic RTSTRUCT import was cancelled.')

          frameDataset = self._datasetForFrame(
            sourceDataset, frameNumber)
          if analysis['contourCountByFrame'].get(str(frameNumber), 0) == 0:
            segmentationNode = self._createEmptyFrameSegmentation(sourceDataset)
          else:
            framePath = os.path.join(
              temporaryDirectory, f'frame_{frameNumber:04d}.dcm')
            pydicom.dcmwrite(
              framePath, frameDataset, write_like_original=False)
            segmentationNode = self._loadTemporaryRTStruct(framePath)
            self._normalizeSegmentIDs(segmentationNode, frameDataset)
            self._ensureAllROIs(segmentationNode, sourceDataset)
          segmentationNode.SetAttribute(
            'dRTImporter.FrameNumber', str(frameNumber))

          if volumeSequence:
            volumeIndexByUID = self._volumeIndexBySOPInstanceUID(volumeSequence)
            referencedUID = analysis['frameToSOPInstanceUID'][str(frameNumber)]
            if volumeIndexByUID is None or referencedUID not in volumeIndexByUID:
              raise RuntimeError(
                f'PET SOP Instance UID for dynamic RT frame {frameNumber} '
                'cannot be mapped to the selected PET sequence.')
            volumeItemNumber = volumeIndexByUID[referencedUID]
            indexValue = volumeSequence.GetNthIndexValue(volumeItemNumber)
          else:
            indexValue = str(frameNumber - 1)
          segmentationSequence.SetDataNodeAtValue(
            segmentationNode, indexValue)
          slicer.mrmlScene.RemoveNode(segmentationNode)

      if volumeBrowser:
        browser = volumeBrowser
      else:
        browser = slicer.mrmlScene.AddNewNodeByClass(
          'vtkMRMLSequenceBrowserNode',
          slicer.mrmlScene.GenerateUniqueName(baseName + ' browser'))

      browser.AddSynchronizedSequenceNode(segmentationSequence)
      browser.SetSaveChanges(segmentationSequence, True)
      browser.SetOverwriteProxyName(segmentationSequence, True)
      browser.SetAttribute(
        'dRTImporter.SequenceNodeID', segmentationSequence.GetID())

      proxySegmentation = browser.GetProxyNode(segmentationSequence)
      if proxySegmentation:
        proxySegmentation.SetAttribute(
          'dRTImporter.DynamicRTStruct', '1')
        proxySegmentation.SetAttribute(
          'dRTImporter.SequenceNodeID', segmentationSequence.GetID())
        proxySegmentation.SetAttribute(
          'dRTImporter.FrameToSOPInstanceUID',
          json.dumps(analysis['frameToSOPInstanceUID']))
        if analysis['referencedSeriesInstanceUID']:
          proxySegmentation.SetAttribute(
            'dRTImporter.ReferencedSeriesInstanceUID',
            analysis['referencedSeriesInstanceUID'])
        if volumeBrowser and volumeSequence:
          volumeProxy = volumeBrowser.GetProxyNode(volumeSequence)
          if volumeProxy:
            proxySegmentation.SetReferenceImageGeometryParameterFromVolumeNode(
              volumeProxy)
        self.addSeriesInSubjectHierarchy(loadable, proxySegmentation)

      if (hasattr(slicer.modules, 'sequences')
          and slicer.modules.sequences.widgetRepresentation()):
        slicer.modules.sequences.widgetRepresentation().setActiveBrowserNode(
          browser)
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

    frameMapJson = sourceSequence.GetAttribute('dRTImporter.FrameToSOPInstanceUID')
    if not frameMapJson:
      raise RuntimeError('The dynamic RTSTRUCT has no frame-to-SOP reference map.')
    frameMap = json.loads(frameMapJson)
    frameNumbers = sorted(int(value) for value in frameMap.keys())
    referencedUIDs = [frameMap[str(frameNumber)] for frameNumber in frameNumbers]

    petUIDs = dRTImporterPluginClass._sequenceFrameSOPInstanceUIDs(petSequence)
    if petUIDs is None:
      raise RuntimeError(
        'The selected PET sequence does not expose one source SOP Instance UID '
        'per temporal frame.')
    if len(petUIDs) != len(referencedUIDs) or set(petUIDs) != set(referencedUIDs):
      raise RuntimeError(
        'The dynamic RTSTRUCT temporal SOP references do not match the '
        'selected PET sequence.')
    petIndexByUID = {uid: index for index, uid in enumerate(petUIDs)}

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

    alignedSequence = slicer.mrmlScene.AddNewNodeByClass(
      'vtkMRMLSequenceNode',
      slicer.mrmlScene.GenerateUniqueName(sourceSequence.GetName() + ' PET-aligned'))
    alignedSequence.CopySequenceIndex(petSequence)
    for attributeName in (
        'dRTImporter.DynamicRTStruct',
        'dRTImporter.Convention',
        'dRTImporter.SourceFile',
        'dRTImporter.FrameToSOPInstanceUID',
        'dRTImporter.ReferencedSeriesInstanceUID'):
      value = sourceSequence.GetAttribute(attributeName)
      if value is not None:
        alignedSequence.SetAttribute(attributeName, value)

    for frameNumber in frameNumbers:
      referencedUID = frameMap[str(frameNumber)]
      petItemIndex = petIndexByUID[referencedUID]
      indexValue = petSequence.GetNthIndexValue(petItemIndex)
      alignedSequence.SetDataNodeAtValue(sourceNodeByFrame[frameNumber], indexValue)

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
    petBrowser.SetOverwriteProxyName(alignedSequence, True)
    petProxy = petBrowser.GetProxyNode(petSequence)
    if petProxy is not None:
      proxyNode.SetReferenceImageGeometryParameterFromVolumeNode(petProxy)
    proxyNode.SetAttribute('dRTImporter.SequenceNodeID', alignedSequence.GetID())
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
