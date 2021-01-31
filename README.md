# fft_three_ways

FFT in 'double', 'float' and fixed-point integer.

This is so I can assess relaative performance and relative error

It generates signed 24-bit random complex values then runs an FFT on them.

## Building

If you have gcc and make installed, then Just type "make" 

## Results

This is for a 16384 point complex FFT.

```
Time taken
==========
Double DFT            0.978966000 secs
Double FFT            0.000853900 secs
Float FFT             0.001034500 secs
Integer FFT           0.001339300 secs
```
