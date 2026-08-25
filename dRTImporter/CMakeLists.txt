#-----------------------------------------------------------------------------
set(MODULE_NAME dRTImporter)

#-----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  dRTImporterLib/__init__.py
  dRTImporterLib/dRTImporterPlugin.py
  dRTImporterLib/dRTExporter.py
)

set(MODULE_PYTHON_RESOURCES
)

set(MODULE_PYTHON_RESOURCES
)

#-----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
)
