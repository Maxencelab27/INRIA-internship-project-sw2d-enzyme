# Aliases
$SW2D_BIN="$SW2D_PATH\build\bin"

function sw2d  # Function to run the Shallow Water 2D suite.
{    
    param(
    [ValidateSet('clean','debug','minsizerel','release','relwithdebinfo')][String]$build='none',
    [ValidateSet('gcc','intel','msvc')][String]$compiler='intel',
    [switch]$doc,
    [switch]$graphics,
    [switch]$help,
    [ValidateSet('converter','default','porosity_mapper','postprocess','solver','zone')][String]$run='none',
    [ValidateSet('all','clean','converter','ddp','dip','liswes','output_converter','porosity_mapper','postprocess','solver','sp','swes','unit','zone')][String]$test='none'
    )

    $path = "."

    # Default config
    if($build -eq 'none' -And !$doc -And !$help -And $run -eq 'none' -And $test -eq 'none')
    {  # The default 'sw2d' runs converter + solver 
        $run = 'default'
    }

    # Sanity check
    if (!(Test-Path $path -PathType Any))
    {
        Write-Output ''
        Write-Output "Error: directory '$path' does not exist!"
        Write-Output ''
        break
    }

    # Build
    switch($build)
    {
        "clean" # Remove 'build' folder
        {
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH
            Remove-Item "build" -Force -Recurse
            Set-Location -Path $curDir
        }
        "debug" # Compile with Debug type
        {
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH
            cmd.exe /c .\recipes\bld_$compiler.bat DEBUG
            Set-Location -Path $curDir
        }
        "minsizerel" # Compile with MinSizeRel type
        {
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH
            cmd.exe /c .\recipes\bld_$compiler.bat MINSIZEREL
            Set-Location -Path $curDir
        }
        "release" # Compile with Release type
        {
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH
            cmd.exe /c .\recipes\bld_$compiler.bat RELEASE
            Set-Location -Path $curDir
        }
        "relwithdebinfo" # Compile with RelWithDebInfo type
        {
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH
            cmd.exe /c .\recipes\bld_$compiler.bat RELWITHDEBINFO
            Set-Location -Path $curDir
        }
        "none" # Do nothing
        {}
        default
        {
            Write-Output ''
            Write-Output 'Unknown SW2D build option : '$build'. Please run *sw2d -help* for more information...'
            Write-Output ''
            break
        }
    }

    # Docs
    if($doc) # Run html doc thanks to doxygen
    {
        $curDir = Get-Location
        Set-Location -Path $SW2D_PATH/doc
        cmd.exe /c doxygen doxyfile
        Invoke-Expression  html/index.html
        Set-Location -Path $curDir
        break
    }

    # Help
    if($help) # Print help related to SW2D 
    {
        Get-Help sw2d -detailed
        break
    }

    # Run
    if (($run -ne 'none') -And (-Not (Test-Path $path\input -PathType Any)) -And ($run -ne 'converter'))
    { # Sanity check
        Write-Output ''
        Write-Output "Error: no input files in path $path"
        Write-Output ''
        break
    }
    else
    {
        switch ($run) {
            "converter" {
                cmd.exe /c $SW2D_BIN\sw2dConverter.exe $path
            }
            "default" { # converter + solver
                cmd.exe /c $SW2D_BIN\sw2dConverter.exe $path
                cmd.exe /c $SW2D_BIN\sw2dSolver.exe $path
            }
            "postprocess" {
                cmd.exe /c $SW2D_BIN\sw2dPostProcess.exe $path
            } 
            "porosity_mapper" {
                cmd.exe /c $SW2D_BIN\sw2dPorosityMapper.exe $path
            }
            "solver" {
                cmd.exe /c $SW2D_BIN\sw2dSolver.exe $path
            }
            "zone" {
                cmd.exe /c $SW2D_BIN\sw2dZone.exe $path
            }
            "none" { # Do nothing
            }
            default {
                Write-Output ''
                Write-Output 'Unknown SW2D executable option '$run'! Please run *sw2d -help* for more information...'
                Write-Output ''
                break
            }
        }
    }

    # Test
    switch ($test) {
        "all" { # Run converter, porosity mapper, solver and zone mr-tests
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH/test/merge/converter/
            cmd.exe /c python run_mr-tests.py
            Write-Output ''
            Read-Host -Prompt "Press any key to continue with postprocess tests"
            Set-Location -Path $SW2D_PATH/test/merge/postprocess/
            cmd.exe /c python run_mr-tests.py
            Write-Output ''
            Read-Host -Prompt "Press any key to continue with porosity mapper tests"
            Set-Location -Path $SW2D_PATH/test/merge/porosity_mapper/
            cmd.exe /c python run_mr-tests.py
            Write-Output ''
            Read-Host -Prompt "Press any key to continue with zone tests"
            Set-Location -Path $SW2D_PATH/test/merge/zone/
            cmd.exe /c python run_mr-tests.py
            Write-Output ''
            Read-Host -Prompt "Press any key to continue with solver tests"
            Set-Location -Path $SW2D_PATH/test/merge/solver/
            cmd.exe /c python run_mr-tests.py
            Set-Location -Path $curDir
        }
        "clean" { # Clean mr-tests
            Remove-Item -Path $SW2D_PATH/test -Force -Recurse
            Write-Output ''
            Write-Output "Cleaning test repository..."
            cmd.exe /c git checkout -- $SW2D_PATH/test
        }
        "converter" { # Run converter mr-tests
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH/test/merge/converter/
            cmd.exe /c python run_mr-tests.py
            Set-Location -Path $curDir
        }
        "ddp" { # Run all DDP mr-tests 
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH/test/merge/solver/ddp/
            Copy-Item -Path "../run_mr-tests.py" -Destination "."
            cmd.exe /c python3 run_mr-tests.py
            Set-Location -Path $curDir
        }
        "dip" { # Run all DIP mr-tests 
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH/test/merge/solver/dip/
            Copy-Item -Path "../run_mr-tests.py" -Destination "."
            cmd.exe /c python run_mr-tests.py
            Set-Location -Path $curDir
        }
        "liswes" { # Run all liswes mr-tests 
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH/test/merge/solver/liswes/
            Copy-Item -Path "../run_mr-tests.py" -Destination "."
            cmd.exe /c python run_mr-tests.py
            Set-Location -Path $curDir
        }
        "output_converter" { # Run solver mr-tests
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH/test/merge/output_converter/
            cmd.exe /c python run_mr-tests.py
            Set-Location -Path $curDir
        }
        "postprocess" { # Run postprocess mr-tests
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH/test/merge/postprocess/
            cmd.exe /c python run_mr-tests.py
            Set-Location -Path $curDir
        }
        "porosity_mapper" { # Run porosity mapper mr-tests
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH/test/merge/porosity_mapper/
            cmd.exe /c python run_mr-tests.py
            Set-Location -Path $curDir
        }
        "solver" { # Run solver mr-tests
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH/test/merge/solver/
            cmd.exe /c python run_mr-tests.py
            Set-Location -Path $curDir
        }
        "sp" { # Run all SP mr-tests 
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH/test/merge/solver/sp/
            Copy-Item -Path "../run_mr-tests.py" -Destination "."
            cmd.exe /c python run_mr-tests.py
            Set-Location -Path $curDir
        }
        "swes" { # Run all swes mr-tests 
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH/test/merge/solver/swes/
            Copy-Item -Path "../run_mr-tests.py" -Destination "."
            cmd.exe /c python run_mr-tests.py
            Set-Location -Path $curDir
        }
        "unit" {# Run all unit-tests
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH/
            cmd.exe /c ./build/bin/sw2dUnitTests.exe 
            Set-Location -Path $curDir
        }
        "zone" { # Run all zone mr-tests 
            $curDir = Get-Location
            Set-Location -Path $SW2D_PATH/test/merge/zone/
            cmd.exe /c python run_mr-tests.py
            Set-Location -Path $curDir
        }
        "none" { # Do nothing
        }
        default {
            Write-Output ''
            Write-Output 'Unknown SW2D test option '$test$'! Please run *sw2d -help* for more information...'
            Write-Output ''
            break
        }
    }

    # Graphics
    if($graphics)
    {
        if($run -eq 'solver' -Or $run -eq 'default')
        {  # Single run: run python script and display png files
            cmd.exe /c python $path\graph_scripts\graph_script.py
            Get-ChildItem $path\graph_scripts\*.png | ForEach-Object {
                cmd.exe /c start $_.fullname
            }
        }
        elseif($test -ne 'none' -And $test -ne 'converter')
        {  # mr-tests run: just open the validation pdf file
            cmd.exe /c start $SW2D_PATH\test\mr_validation\validation_doc\ShortValidationDoc.pdf
        }
    }

    <#
        .SYNOPSIS
        Runs the Shallow Water 2D suite.

        .DESCRIPTION
        Runs the Shallow Water 2D suite.
        Takes several types of options:
        sw2d [-b <build_option> ] [-c <compiler_option> ] [-d] [-g] [-h] [-r <run_option>] [-t <test_option>]

        .PARAMETER build
        -b | -build
        Builds all SW2D libraries and binaries.
        <String> = clean           -> remove 'build' folder
        <String> = debug           -> compile with Debug build type
        <String> = minsizerel      -> compile with MinSizeRel build type
        <String> = release         -> compile with Release build type
        <String> = relwithdebinfo  -> compile with RelWithDebInfo build type

        .PARAMETER compiler
        -c | -compiler
        Specifies the compiler when parameter build is used. If the option is not defined, default compiler is Intel.
        <String> = gcc     -> use gcc compiler
        <String> = intel   -> use intel compiler
        <String> = msvc    -> use msvc compiler

        .PARAMETER doc
        -d | -doc
        Run html (doxygen) documentation.

        .PARAMETER graphics
        -g | -graphics
        Displays graphics (or not) after running code.
        If the option is used with a solver single run, displays the simulation graph scripts after the run.
        It the option is used with solver mr-tests, displays the validation doc pdf after running the tests.

        .PARAMETER help
        -h | -help
        Shows this information.

        .PARAMETER run
        -r | -run
        Specifies the SW2D program to be run.
        <String> = converter | solver | postprocess | porosity_mapper | zone -> runs the corresponding binary
        <String> = default                                                   -> runs converter, then solver

        .PARAMETER test
        -t | -test
        Runs some sets of tests.
        <String> = all             -> runs all converter + porosity_mapper + postprocess + solver tests
        <String> = clean           -> remove 'test' folder and restore it from gitlab
        <String> = converter       -> runs all converter mr-tests 
        <String> = output_converter -> runs all output converter mr-tests
        <String> = porosity_mapper -> runs all porosity mapper mr-tests
        <String> = postprocess     -> runs all postprocess mr-tests
        <String> = solver          -> runs all solver mr-tests
        <String> = unit            -> runs all unit tests
        <String> = zone             -> runs all zone mr-tests
        <String> = $tag            -> runs all solver mr-tests with corresponding to $tag in {ddp, dip, liswes, mp, sp, swes}

        .EXAMPLE
         COMMANDS           | WHAT IT DOES
    ------------------------|-----------------------------------------------------------------
    sw2d -b debug           | Builds SW2D targets with Debug build type and Intel compiler
    sw2d -b release -c gcc  | Builds SW2D targets with Release build type and gcc compiler
    sw2d -d                 | Runs doxygen (html) documentation in the sw2d/doc directory
    sw2d -h                 | Displays this help menu
    sw2d -t solver -g       | Runs all the solver mr-tests and display the validation pdf file
    sw2d -run solver        | Runs solver from the current directory

        .LINK
        See SW2D web site: https://sw2d.inria.fr
    #>
}