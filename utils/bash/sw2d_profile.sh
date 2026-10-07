#!/usr/bin/env bash

# SW2D profile (mainly used to define the sw2d function)
## This file should be loaded from your .bash_profile
## It can be used with the sw2d autocompletion tools

# Add binaries in PATH
export PATH="$PATH$SW2D_PATH/build/bin/:$PATH" 
alias plotmesh="python $SW2D_PATH/utils/python/plot/mesh_plot/mesh_all.py"
alias plotmesh_bnd="python $SW2D_PATH/utils/python/plot/mesh_plot/mesh_boundary.py"
## WARNING: should use recent version of getopt function in bash. 
## For Mac Users, use brew install getopt and adapt $PATH

# Build the SW2D function
unset -f sw2d
function sw2d
{
    local build='none'
    local currentpath=$(pwd)
    local doc='none'
    local help='none'
    local run='default'
    local run_sw2d=0

    options=$(getopt -q -o bdhgr:t: --long build,doc,help,graphics,run:, t: -- "$@")
    if [[ $? != 0 ]] ; then echo -e "" ; echo "Unkown SW2D sequence: 'sw2d $@'. Please run 'sw2d --help' for more info." ; echo -e "";
    else
    {
        while true; do
            case "$1" in
                -b|--build) # Build SW2D targets
                    if [ -z "$2" ] ; then build="same"; else build="$2";
                    fi
                    cd $SW2D_PATH
                    if [ ! -d "build" ]; then echo -e "Creating build directory..."; mkdir build;
                    fi
                    cd build
                    case "${build}" in
                        clean)
                            echo -e ""
                            echo -e "Cleaning SW2D targets..."
                            echo -e ""
                            cd ..
                            rm -rf build/
                            mkdir build
                            ;;
                        Coverage)
                            echo -e ""
                            echo -e "Building SW2D targets with Coverage build type..."
                            echo -e ""
                            rm CMakeCache.txt Makefile
                            rm -rf test
                            cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_COVERAGE=ON ..
                            make -j
                            ;;
                        Debug|Minisizerel|Release|RelWithDebInfo)
                            echo -e ""
                            echo -e "Building SW2D targets with ${build^^} build type..."
                            echo -e ""
                            rm CMakeCache.txt Makefile
                            cmake -DCMAKE_BUILD_TYPE=${build^^} ..
                            make -j
                            ;;
                        same)
                            echo -e ""
                            echo -e "Building SW2D targets with same build type..."
                            echo -e ""
                            make -j
                            ;;
                        *)
                            echo -e ""
                            echo -e "Unknown '$build' option for target '$1'. Please run 'sw2d --help' for more info."
                            echo -e ""
                            break;;
                    esac
                    cd $currentpath
                    break;;
                -d|--doc) # Doc option (to be used through alias): read examples catalog or validation pdf
                    doc="$2"
                    if [ -z "$2" ] ; then doc="all"; else doc="$2";
                    fi
                    cd $SW2D_PATH/doc
                    case "${doc}" in
                        all)
                            doxygen doxyfile
                            open html/index.html
                            cd pdf/src/
                            pdflatex --interaction=batchmode developer-main.tex
                            pdflatex --interaction=batchmode developer-main.tex
                            pdflatex --interaction=batchmode user-main.tex
                            pdflatex --interaction=batchmode user-main.tex
                            mv developer-main.pdf ../sw2d_developer.pdf
                            mv user-main.pdf ../sw2d_user.pdf
                            rm -rf _markdown* *.aux *.lua *.out *.log *.toc
                            open ../sw2d_developer.pdf
                            open ../sw2d_user.pdf
                            ;;
                        clean)
                            cd ..
                            rm -rf doc
                            git restore doc
                            ;;
                        html)
                            doxygen doxyfile_with_version
                            open html/index.html
                            ;;
                        pdf)
                            cd pdf/src/
                            pdflatex --interaction=batchmode developer-main.tex
                            pdflatex --interaction=batchmode developer-main.tex
                            pdflatex --interaction=batchmode user-main.tex
                            pdflatex --interaction=batchmode user-main.tex
                            mv developer-main.pdf ../sw2d_developer.pdf
                            mv user-main.pdf ../sw2d_user.pdf
                            rm -rf *_markdown *.aux *.lua *.out *.log *.toc
                            open ../sw2d_developer.pdf
                            open ../sw2d_user.pdf
                            ;;
                        *)
                            echo -e ""
                            echo -e "Unknown '$doc' option for target '$1'. Please run 'sw2d --help' for more info."
                            echo -e ""
                            break;;
                    esac
                    cd $currentpath
                    break;;
                -h|--help) # Help
                    echo -e ""
                    case "${help}" in
                        none)
                            sw2d_help
                            ;;
                    esac
                    echo -e ""
                    break;;
                -g|--graphics) # Display output plots
                    echo -e ""
                    python graph_scripts/graph_script.py
                    open graph_scripts/*.png
                    echo -e ""
                    break;;
                -r|--run) # Execution option: defines the type of SW2D binary to run.
                    run="$2"
                    case "${run}" in
                        converter|porosity_mapper|postprocess|solver|zone)
                            run_sw2d=1
                            ;;
                        *)
                            echo -e ""
                            echo -e "Unknown '$run' option for target '$1'. Please run 'sw2d --help' for more info."
                            echo -e ""
                            break
                            ;;
                    esac
                    shift 2;;
                -t|--test) # Testing mr_tests
                    if [ -z "$2" ] ; then test="here"; else test="$2";
                    fi
                    case "${test}" in
                        clean)
                            rm -rf $SW2D_PATH/test
                            echo -e ""
                            echo -e "Cleaning test repository..."
                            git checkout -- $SW2D_PATH/test
                            cd $currentpath
                            ;;
                        converter|porositymapper|postprocess|zone)
                            echo -e ""
                            cd $SW2D_PATH/test/merge/$test/
                            python test$test.py
                            cd $currentpath
                            open $SW2D_PATH/doc/pdf/sw2d_mergeRequestDoc.pdf
                            ;;
                        full)
                            echo -e ""
                            cd $SW2D_PATH/test/merge/
                            python buildValidationDoc.py
                            cd $currentpath
                            open $SW2D_PATH/doc/pdf/sw2d_validationDoc.pdf
                            ;;
                        solver)
                            echo -e ""
                            cd $SW2D_PATH/test/merge/solver/
                            python testSolver.py
                            cd $currentpath
                            open $SW2D_PATH/doc/pdf/sw2d_mergeRequestDoc.pdf
                            ;;
                        ddp|dip|liswes|sp|swes)
                            echo -e ""
                            cd $SW2D_PATH/test/merge/solver/_in_progress/$test
                            cp $SW2D_PATH/test/merge/solver/testSolver.py .
                            python $SW2D_PATH/test/merge/solver/testSolver.py .
                            rm testSolver.py
                            cd $currentpath
                            open $SW2D_PATH/doc/pdf/sw2d_mergeRequestDoc.pdf
                            ;;
                        unit)
                            echo -e ""
                            cd $SW2D_PATH
                            ./build/bin/sw2dUnitTests
                            #gcovr -r src/common . --txt -> build/coverage.txt
                            gcovr -r src . --filter src/common --filter src/sw2dSolver/common --txt > build/coverage.txt
                            open build/coverage.txt
                            cd $currentpath
                            ;;
                        here)
                            echo -e ""
                            cp $SW2D_PATH/test/merge/solver/testSolver.py .
                            python testSolver.py
                            rm testSolver.py
                            ;;
                        *)
                            echo -e ""
                            echo -e "Unknown '$test' option for target '$1'. Please run 'sw2d --help' for more info."
                            echo -e ""
                            break;;
                    esac
                    case "${3}" in
                        '')
                            break;;
                    esac
                    shift 2;;
                --) 
                    shift
                    break;;
                "") # Default config
                    run_sw2d=1
                    break;;
                * | ?)
                    echo -e ""
                    echo -e "Unknown sw2d target: '$1'. Please run 'sw2d --help' for more info."
                    echo -e ""
                    break;;
            esac
        done

        if [[ $run_sw2d == 1 ]]
        then
        {
            if [ -d input ] ; then
            {
                if [[ $run == 'converter' ]]
                then 
                    $SW2D_PATH/build/bin/sw2dConverter .
                elif [[ $run == 'porosity_mapper' ]]
                then 
                    $SW2D_PATH/build/bin/sw2dPorosityMapper .
                elif [[ $run == 'postprocess' ]]
                then 
                    $SW2D_PATH/build/bin/sw2dPostProcess .
                elif [[ $run == 'solver' ]]
                then 
                    $SW2D_PATH/build/bin/sw2dSolver .
                elif [[ $run == 'zone' ]]
                then 
                    $SW2D_PATH/build/bin/sw2dZone .
                elif [[ $run == 'default' ]]
                then
                    $SW2D_PATH/build/bin/sw2dConverter .
                    $SW2D_PATH/build/bin/sw2dSolver .
                else # Sanity check
                {
                    echo -e ""
                    echo -e "Unknown '$run' option for target '$1'. Please run 'sw2d --help' for more info."
                    echo -e ""
                }
                fi
            }
            else # Sanity check
            {
                echo -e ""
                echo -e "No input folder here. Please run 'sw2d --help' for more info."
                echo -e ""
            }
            fi
        }
        fi
    }
    fi
}

# Autocompletion
source $SW2D_PATH/utils/bash/bash-completion.d/sw2d.sh

# Help function
sw2d_help(){
>&2 cat << EOF
Help: sw2d

sw2d [-b [<build_options>]] [-d [<doc_options>]] [-h] [-g] [-r <run_options>] [-t <test_options>]

1. sw2d options
        These options should be run after the *sw2d* command only (no alias). Short and long options are available.
        -b | --build     optional input in {clean, Debug, Minisizerel, RelWithDebInfo, Release}. 
                           Builds SW2D targets (from anywhere). With no input, cmake files are conserved.
        -d | --doc       optional input in targets {all, html, pdf}. Runs doxygen (html) or latex (pdf) documentation. Default value is all.
        -h | --help      displays this help
        -g | --graphics  runs and displays Python graph scripts (in local project folder only)
        -r | --run       mandatory input in {converter, porosity_mapper, postprocess, solver, zone}
        -t | --test      mandatory input in {clean, converter, ddp, dip, full, liswes, porosity_mapper, postprocess, solver, sp, swes, unit, zone}. 
                           Runs (part of) the merge tests thanks to the python script.

2. Examples
    COMMANDS              | WHAT IT DOES
    ----------------------|--------------------------------------------------------------------------------------------------------
    sw2d --build <type>   | Builds SW2D targets from any directory (and brings you back in) with optional <type> cmake built type
    sw2d -d               | Runs doxygen (html) documentation in the sw2d/doc directory
    sw2d --help           | Displays this help menu
    sw2d -g               | Runs graph_script/graph_script.py and displays corresponding png figures.
    sw2d --run converter  | Runs the converter in the current directory
    sw2d --test swes      | Runs all merging tests in the swes directory
EOF
}