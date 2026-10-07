### Using SW2D shortcuts under bash (MacOS, Linux)

These utilities and shortcuts are provided to make using the code easier, but their functionality is not guaranteed. The authors do not guarantee that they will be able to help users use these tools.


#### The sw2d function 

The file `utils/bash/sw2d_profile.sh` should be sourced in your bash terminal. This can be done thanks to the addition of the following lines in your `.bash_profile`: <br> <br> 

`export SW2D_PATH=<PATH-TO-SW2D>`
`if [ -f $SW2D_PATH/utils/bash/sw2d_profile.sh ]; then
  export PATH="/opt/homebrew/opt/gnu-getopt/bin:$PATH"
  . $SW2D_PATH/utils/bash/sw2d_profile.sh`
`fi`

Note that gnu-getopt should be installed (with brew for MacOS users), in order to manage long option names in sw2d call. Make sure that the path is compatible in your bash profile:

`export PATH="/opt/homebrew/opt/gnu-getopt/bin:$PATH"`

Don't forget to close and open your terminal console after saving your profile file.

The environment variable `SW2D_PATH` is used in many places, please make sure the path is correctly indicated.
The `sw2d_profile.sh` file creates the `sw2d` command that can be run with options and directory arguments. See `sw2d -h` for more information.

#### sw2d arguments

Run `sw2d -h` for more information.

#### Autocompletion 

The file `utils/bash/bash-completion.d/sw2d.sh` is used for autocompletion, including dedicated arguments. It is meant to be sourced from the `sw2d_profile` file, or copied in your `etc/bash_completion.d/` folder.
With brew, the default bash completion adress is `/opt/homebrew/etc/bash_completion.d/`. **This feature is a WIP**, with no warranty. In particular, it should use bash version >=4.