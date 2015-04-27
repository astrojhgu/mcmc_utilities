#!/usr/bin/env python
import scipy

x=-10
k=4
b=3
while x<100:
    y=x*k+b
    y_obs=scipy.random.normal(0,10.89)+y
    print(x,y_obs)
    x+=.1
