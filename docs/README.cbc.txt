The solver "cbc" uses Cbc (Coin-or branch and cut; see
https://github.com/coin-or/Cbc) to solve mixed-integer linear programming
problems.

To use Cbc from AMPL, put the solver executable in a directory in your
usual PATH and run:

     option solver cbc;
     solve;

AMPL invokes "cbc stub -AMPL", where stub.nl contains the model, and
reads the resulting stub.sol file automatically.

Control solver options with AMPL's cbc_options option or the cbc_options
environment variable. Run "cbc -=" to see the available options.

INSTALLING
==========

On Windows, put cbc.exe in a directory in your PATH.
On Linux and macOS, keep cbc executable (chmod +x cbc) and put it in a
directory in your PATH.

For support, contact support@ampl.com.
