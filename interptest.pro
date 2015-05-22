; Written May 21 for SURP 2015 project
; Tests vectorized vs looped interpolation

; Create array that serves as values to interpolate
x = findgen(1e4)
; Create array at which we wish to know the values
y = make_array(1e5,/index,start=0, increment = 0.1)

; Time vectorized interpolation and print results
timer, /start, number = 1
R1 = interpol(x,x,y,/spline)
timer, /stop, number = 1
timer, /print, number = 1

; Time looped interpolation and print results
timer, /start, number = 2
R2 = fltarr(leng(y))
for i=0, leng(y)-1 do begin
R2[i] = interpol(x,x,y[i],/spline)
endfor
timer, /stop, number = 2
timer, /print, number = 2

end
