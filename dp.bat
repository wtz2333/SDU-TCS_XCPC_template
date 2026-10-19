:loop
gen.exe
std.exe < data.in > std.out
brute.exe < data.in > brute.ans
fc std.out brute.ans
if not errorlevel 1 goto loop
pause