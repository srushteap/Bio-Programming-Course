#!/usr/bin/env python3
# Name: Srushti Patil(spatil5)
# Group Members: none

"""
Build a Person object and have it introduce itself.

input: a string of arbitrary length, which is used to name the new person object
output: greeting printed to screen
"""

class Person:
    def __init__(self,name,pet): #creates objects for Person class
        self.myName = name
        self.myPet  = pet
        
    def introduce (self): #introduce method prints an introduction using objects from Person class
        print ("Hi there, I am {0}, and I like {1}s".format(self.myName,self.myPet))


# put your new code here.
print("What is your name?: ") #asks user for name and favorite kind of pet
name = input()
print("What is your favorite kind of pet? (singular): ")
pet = input()
newPerson = Person(name, pet) #creates a new name of class Person with the inputs as the objects
newPerson.introduce() #sends introduce() method message to print out introduction using inputs