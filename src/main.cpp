#include "command_line.hpp"
#include "fusion.hpp"
#include "progress.hpp"
#include "pwmvs.hpp"
#include "pwmvs_controller.hpp"
#include "workspace_io_eth3d.hpp"
#include "workspace_io_openmvg.hpp"

#include <openMVG/sfm/sfm.hpp>

int main( int argc, char** argv )
{
  // Handle command line
  bool            geometric;
  std::string     input_file_path, output_folder_path;
  PWMVS::Options  pwmvs_options;
  Fusion::Options fusion_options;
  handle_global_parameters( input_file_path, output_folder_path, geometric, argc, argv );
  handle_pwmvs_arguments( pwmvs_options, argc, argv );
  handle_fusion_arguments( fusion_options, argc, argv );

  // Initialize workspace
  std::shared_ptr<Workspace> workspace = std::make_shared<Workspace>();
  InitializeWorkspaceOpenMVG( input_file_path, output_folder_path, *workspace, 0.5 );
  workspace->work_path = workspace->root_path + "/pwmvs";
  workspace->initialize();
  std::cout << "Initialize done" << std::endl;

  // Perform pwmvs computation
  Controller<PWMVS> pwmvs_controller( workspace );
  pwmvs_controller.pwmvs_options = pwmvs_options;
  ConsoleProgress pwmvs_progress;
  if ( !pwmvs_controller.run( geometric, &pwmvs_progress ) )
    return EXIT_FAILURE;

  // Perform fusion
  Fusion fusion( workspace );
  fusion.options = fusion_options;
  fusion.run( geometric, &pwmvs_progress );

  // Export points
  std::cout << "Exporting points" << std::endl;
  ExportPoints( fusion.points, fusion.normals, fusion.colors, stlplus::create_filespec( workspace->work_path, "fused", "ply" ) );

  return EXIT_SUCCESS;
}