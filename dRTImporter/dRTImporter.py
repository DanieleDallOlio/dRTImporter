import slicer


class dRTImporter:
  """Hidden scripted module that registers the dynamic RTSTRUCT DICOM plugin."""

  def __init__(self, parent):
    parent.title = "DICOM Dynamic RTSTRUCT Import Plugin"
    parent.categories = ["Developer Tools.DICOM Plugins"]
    parent.dependencies = ["DICOM", "Sequences", "Segmentations", "DicomRtImportExport"]
    parent.contributors = [
      "Daniele Dall'Olio (University of Bologna)"
    ]
    parent.helpText = (
      "DICOM plugin that recognizes dynamic RT Structure Sets and "
      "loads them as segmentation sequences. Static RT Structure Sets are "
      "intentionally ignored."
    )
    parent.acknowledgementText = (
      "Developed for dynamic PET research workflows in 3D Slicer. "
      "RT Structure Set conversion is provided by SlicerRT."
    )
    parent.hidden = True

    from dRTImporterPlugin import dRTImporterPluginClass

    try:
      slicer.modules.dicomPlugins
    except AttributeError:
      slicer.modules.dicomPlugins = {}

    slicer.modules.dicomPlugins['dRTImporterPlugin'] = dRTImporterPluginClass
