#!/usr/bin/env python3
# Name: Srushti Patil(spatil5)
# Group Members: none

"""
Build a Person object and have it introduce itself.

input: a string of arbitrary length, which is used to name the new person object
output: greeting printed to screen
"""

class Person:
    """ Build a new Person objects with all of the attributes given at instantiation time. """
    def __init__(self,name, username, major, studentType, reason, bioInterest, experience): #creates objects for Person class
        """ Save all of the instance arguments. """
        self.myName = name
        self.myMajor = major
        self.myUsername = username
        self.myStudentType = studentType
        self.myReason = reason
        self.myBioInterest = bioInterest
        self.myExperience = experience

    def introduce (self): #prints the introduction using object atributes
        ''' Print the introduction using all of the object atributes '''
        print ("Hi there, I am {0}.".format(self.myName))
        print ("My username is {0}.".format(self.myUsername))
        print ("I am a {0}.".format(self.myStudentType))
        print ("My major is {0}.".format(self.myMajor))
        print (self.myReason)
        print (self.myBioInterest)
        print (self.myExperience)

name = 'Srushti Patil' #new variables to describe user
username = 'spatil5'
studentType = 'student'
major = 'Biomolecular Engineering'
reason = 'I\'m taking this class because it is a major requirement. I\'m excited to learn how to apply coding in a bio research setting.'
bioInterest = 'I am interested in stem cells.'
experience = 'I took a few coding classes in highschool and I also took CSE 20 here at UCSC so I am familiar with java and python.'
srushti = Person(name, username, major, studentType, reason, bioInterest, experience) #creates a new name of class Person with the above variables as objects
srushti.introduce() #sends message to introduce() to write introduction using above variables
