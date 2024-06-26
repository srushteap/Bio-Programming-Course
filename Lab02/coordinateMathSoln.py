#!/usr/bin/env python3 
# Name: Srushti Patil(spatil5)
# Group Members: Andres, Dylan, Valerie

'''
Calculates the bond length and angles from three sets of user inputted atomic coordinates. 

Example:
 Input: C = (39.447, 94.657, 11.824) N = (39.292, 95.716, 11.027) Ca = (39.462, 97.101, 11.465)
Output: N-C bond length = 1.33
        N-Ca bond length = 1.46
        C-N-Ca bond angle = 124.0
'''

import math
class Triad:
    """
    Calculate angles and distances among a triad of points.
 
    Author: David Bernick
    Date: March 21, 2013
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()
 
    Required Modules: math
    initialized: 3 positional tuples representing Points in n-space
             p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,1,0) )
    attributes: p,q,r the 3 tuples representing points in N-space
    methods:  angleP(), angleR(), angleQ() angles measured in radians
          dPQ(), dPR(), dQR() distances in the same units of p,q,r
 
    """
 
    def __init__(self,p,q,r) :
        """ Construct a Triad. 
        
        Example construction:
            p1 = Triad( p=(1.,0.,0.), q=(0.,0.,0.), r=(0.,0.,0.) ). 
        """
        self.p = p
        self.q = q
        self.r = r
        
# private helper methods
    def d2 (self,a,b) : # calculate squared distance of point a to b
        return float(sum((ia-ib)*(ia-ib)  for  ia,ib in zip (a,b)))
    
    def dot (self,a,b) : # dotProd of standard vectors a,b
        return float(sum(ia*ib for ia,ib in zip(a,b)))
    
    def ndot (self,a,b,c) : # dotProd of vec. a,c standardized to b
        return float(sum((ia-ib)*(ic-ib) for ia,ib,ic in zip (a,b,c)))
    
# calculate lengths(distances) of segments PQ, PR and QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return math.sqrt(self.d2(self.p,self.q))
    
    def dPR (self):
        """ Provides the distance between point p and point r """
        return math.sqrt(self.d2(self.p,self.r))
    
    def dQR (self):
        """ Provides the distance between point q and point r """
        return math.sqrt(self.d2(self.q,self.r))
    
    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return math.acos(self.ndot(self.q,self.p,self.r) / 
                         math.sqrt(self.d2(self.q,self.p) * self.d2(self.r,self.p)))
    
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return math.acos(self.ndot(self.p,self.q,self.r) /
                         math.sqrt(self.d2(self.p,self.q) * self.d2(self.r,self.q)))
 
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return math.acos(self.ndot(self.p,self.r,self.q) /
                         math.sqrt(self.d2(self.p,self.r) * self.d2(self.q,self.r)))

def main():
    '''Asks user to input coordinates of atom and prints bond lengths and angle.'''
    
    '''input string'''
    coord = input("Atomic coordinates? (one line only in format: atom = (x,y,z) )")
    splCoord = coord.split(")")
    splCoord.pop()
    dict_coords = {}
    
    '''splits input string and creates dictionary'''
    for i in splCoord:
        spl2 = i.split('= (')
        dict_coords[spl2[0].strip()] = spl2[1].strip()
     
    '''dictionary values to float tuple'''
    dict_coords = {key: tuple(map(float, value.split(','))) for key, value in dict_coords.items()}

    '''creates new name of class Triad with user input values'''
    newTriad = Triad(dict_coords['C'], dict_coords['N'], dict_coords['Ca'])
    
    '''calculates bond length by calling bond length method for corresponding coordinate and prints'''
    lengthNC = newTriad.dPQ() 
    print(f'N-C bond length = {lengthNC:.2f}')
    lengthNCa = newTriad.dPR() 
    print(f'N-Ca bond length = {lengthNCa:.2f}')
    
    '''calculates bond angle by calling bond angle method for corresponding coordinate and prints'''
    angleCNCa = newTriad.angleQ() * 180 / math.pi #converts angle from radians to degrees
    print(f'C-N-C bond angle = {angleCNCa:.1f}')
    

main()