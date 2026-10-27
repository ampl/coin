The solver "ipopt" uses Ipopt (Interior Point OPTimizer; see
https://github.com/coin-or/Ipopt) to solve continuous nonlinear
programming problems.
Normally ipopt is invoked by AMPL's solve command, which gives the
invocation

     ipopt stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
ipopt writes a stub.sol file for use by ampl's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver ipopt;
     solve;

You can control ipopt by setting the environment variable ipopt_options
appropriately (either by using ampl's option command, or by using the
shell's set and export commands before you invoke ampl).  You can put
one or more (white-space separated) phrases in $ipopt_options.  To see
the possibilities, invoke

        ipopt -=

This build solves its linear systems with MUMPS.

INSTALLING
==========

On MS Windows systems, ipopt.exe and libipoptfort.dll must appear
somewhere in your usual search $PATH (or in the current directory).

On Linux and MacOSX systems, mark the binary executable
(chmod +x ipopt) and put it in a directory in your usual $PATH
(or in the current directory).

If you have questions about or find bugs with this stuff,
please contact:

     AMPL Support
     support@ampl.com
