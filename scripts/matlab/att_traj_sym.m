syms t real
syms amp ampPercent freq freqOffset ampOffset

yd = amp * (1 - ampPercent + ampPercent * cos(freq * t + freqOffset)*cos(freq * t + freqOffset)) + ampOffset
dyd = diff(yd,t)
d2yd = diff(dyd,t)

