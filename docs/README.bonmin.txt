The solver "bonmin" uses Bonmin (Basic Open-source Nonlinear Mixed
INteger programming; see https://github.com/coin-or/Bonmin), together
with Cbc and Ipopt, to solve mixed-integer nonlinear programming
problems; it is exact for convex problems and a heuristic otherwise.
Normally bonmin is invoked by AMPL's solve command, which gives the
invocation

     bonmin stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
bonmin writes a stub.sol file for use by ampl's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver bonmin;
     solve;

You can control bonmin by setting the environment variable bonmin_options
appropriately (either by using ampl's option command, or by using the
shell's set and export commands before you invoke ampl).  You can put
one or more (white-space separated) phrases in $bonmin_options.  To see
the possibilities, invoke

        bonmin -=

INSTALLING
==========

On MS Windows systems, bonmin.exe and libipoptfort.dll must appear
somewhere in your usual search $PATH (or in the current directory).

On Linux and MacOSX systems, mark the binary executable
(chmod +x bonmin) and put it in a directory in your usual $PATH
(or in the current directory).

If you have questions about or find bugs with this stuff,
please contact:

     AMPL Support
     support@ampl.com
