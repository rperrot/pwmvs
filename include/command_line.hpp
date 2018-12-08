#ifndef _PWMVS_COMMAND_LINE_HPP_
#define _PWMVS_COMMAND_LINE_HPP_

#include "fusion.hpp"
#include "pwmvs.hpp"

#include <string>

/**
 * @brief Display usage of the program then exit
 * 
 * @param pgm Program name 
 */
void usage( const std::string& pgm );

/**
 * @brief Get command line arguments for global parameters 
 * 
 * @param input_path      Input SfM file path 
 * @param output_path     Output pwmvs folder path 
 * @param geometric       Use geometric 
 * @param argc            main arguments count 
 * @param argv            main arguments values 
 */
void handle_global_parameters( std::string& input_path, std::string& output_path, bool& geometric, int argc, char** argv );

/**
 * @brief Get command line arguments for PWMVS 
 * 
 * @param pwmvs_options   Options to set
 * @param argc            main arguments count 
 * @param argv            main arguments values 
 */
void handle_pwmvs_arguments( PWMVS::Options& pwmvs_options, int argc, char** argv );

/**
 * @brief Get command line arguments for fusion 
 * 
 * @param fusion_options  Options to set 
 * @param argc            main arguments count 
 * @param argv            main arguments values 
 */
void handle_fusion_arguments( Fusion::Options& fusion_options, int argc, char** argv );

#endif