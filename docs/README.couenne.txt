The solver "couenne" uses Couenne (Convex Over and Under ENvelopes
for Nonlinear Estimation; see https://github.com/coin-or/Couenne),
together with Bonmin, Cbc and Ipopt, to find global optima of
(nonconvex) mixed-integer nonlinear programming problems.
Normally couenne is invoked by AMPL's solve command, which gives the
invocation

     couenne stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
couenne writes a stub.sol file for use by ampl's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver couenne;
     solve;

You can control couenne by setting the environment variable
couenne_options appropriately (either by using ampl's option command,
or by using the shell's set and export commands before you invoke
ampl).  You can put one or more (white-space separated) phrases in
$couenne_options.  To see the possibilities, invoke

        couenne -=

INSTALLING
==========

On MS Windows systems, couenne.exe and libipoptfort.dll must appear
somewhere in your usual search $PATH (or in the current directory).

On Linux and MacOSX systems, mark the binary executable
(chmod +x couenne) and put it in a directory in your usual $PATH
(or in the current directory).

If you have questions about or find bugs with this stuff,
please contact:

     AMPL Support
     support@ampl.com
