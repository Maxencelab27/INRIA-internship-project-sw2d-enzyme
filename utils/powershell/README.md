### Using SW2D shortcuts under Powershell (Windows)

These utilities and shortcuts are provided to make using the code easier, but their functionality is not guaranteed. The authors do not guarantee that they will be able to help users use these tools.


#### The sw2d function 

The file `utils\powershell\sw2d_profile.ps1` should be sourced in your powershell environment. This can be done thanks to slight modifications of your Powershell profile.
Please check that this file - which adress is given by `echo $PROFILE` - does exist (the adress is generally 

- `$home\Documents\WindowsPowerShell\Microsoft.PowerShell_profile.ps1` 
- or `$home\Documents\PowerShell\Microsoft.PowerShell_profile.ps1`).

If this file does not exist, you may create an empty one with the following:

- `New-Item -Path $PROFILE -ItemType File -Force`

If the profile file already exists, please add the following lines to load the sw2d commands (don't forget to change the path to your SW2D repository):

- `$SW2D_PATH = "<PATH-TO-SW2D>"`
- ` . $SW2D_PATH\utils\powershell\sw2d_profile.ps1`

Don't forget to close and open your powershell console after saving your profile file.

The environment variable `SW2D_PATH` is used in many places, please make sure the path is correctly indicated.
The `sw2d_profile.ps1` file creates the `sw2d` function that can be run with optional and directory arguments. See `sw2d -h` for more information.

#### sw2d arguments

See `sw2d -h` for more information.

#### Autocompletion 

The autocompletion is activated for all sw2d arguments.