# Usage
VaSeBuilder offers users to run the program in different ways, known as runmodes.
VaSeBuilder can be run via "python vase.py" and setting the required general parameters as well as the output mode specific parameters.
VaSeBuilder can also be build a set of FastQ files with spiked in donor reads and variants 

## Output modes
VaSeBuilder can be run in different ways called output modes which can be set via the output-mode parameter. Both the output, as well as the work VaSeBuilder does depends on the selected output mode. Some modes do output a validation set, whereas others do not.

### A-mode


### D-mode

Note that D-mode will be implemented in the future.


### P-mode


### V-mode
V-mode only establishes variant contexts and outputs a variant context file. This mode therefore does not exchange acceptor reads with donor reads. This mode can be helpful for example when you want to inspect established variant contexts when using different input data or options.


## Example commands

VaSeBuilder A-mode
<command>

VaSeBuilder P-mode
<command>

VaSeBuilder V-mode
<command>
