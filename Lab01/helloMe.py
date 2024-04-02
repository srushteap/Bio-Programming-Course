#!/usr/bin/env python3
# Name: Srushti Patil(spatil5)
# Group Members: none

"""
Print out the most common statement in computer programming, Hello World.
""" 
print ("Hello World")
''' this is a single line docstring '''
class Announcer (str): #Announcer is now defined as a kind of python string.
    def printMe (self): #printMe method allows objects of class anouncer to print themselves
        print (self)

student = Announcer ('Hello Srushti') #creates a new name of class announcer. parameter is given to class announcer as its data
dog = Announcer ('bark bark') 
student.printMe() #sends printMe() message to object which then prints out its data
dog.printMe()