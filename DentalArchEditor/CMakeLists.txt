set(MODULE_NAME DentalArchEditor)

# Keep the pure numerical/session helpers beside the installed Python module.
file(GLOB MODULE_LIB_SCRIPTS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "DentalArchEditorLib/*.py")
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  ${MODULE_LIB_SCRIPTS}
)
file(GLOB_RECURSE MODULE_PYTHON_RESOURCES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/*")

slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
)
