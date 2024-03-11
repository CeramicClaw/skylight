# skylight
Plotting celestial objects in the sky

I'm using the following papers for calculating the position of celestial objects:
## Moon
Reda, I. *Solar Eclipse Monitoring for Solar Energy Applications Using the Solar and Moon Position Algorithms*; NREL/TP-3B0-47681; National Renewable Energy Laboratory: Golden, Colorado, March 2010. (http://www.nrel.gov/docs/fy10osti/47681.pdf)

### Errata in Moon Position Algorithm
The referenced paper above contains several frustrating mistakes:

1. Equation 26: The equation as written ($a sin$) is $arcsin$
2. Equation 49: The denominator of the equation ($cos \delta - y*sin \pi * cos H$) incorrectly uses $y$ when it should be using $x$. The denominator should be: $cos \delta - x*sin \pi * cos H$

## Sun
Reda, I; Andreas, A. *Solar Position Algorithm for Solar Radiation Applications*; NREL/TP-560-34302; National Renewable Energy Laboratory: Golden, Colorado, January 2008. (https://www.nrel.gov/docs/fy08osti/34302.pdf)