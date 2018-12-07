# -*- coding: utf-8 -*-
'''
Created on Mon Oct 19 17:30:27 2015
Provides some help for nipype workflows 
(to be used within connect statements mostly)
@author: Pierre-Yves Herve
'''


def makeFStringElementFromFnameList(inlist, replacefrom, replaceto, getBasename=True):
    '''
    Create a list of format string or string basd on a list of 
    file names. (To use with nipype Rename or provide output fname/postfix 
    that depends on the input fname patterns.)
    '''
    import os.path as op
    
    if getBasename:
        return [op.basename(s.replace(replacefrom, replaceto)) for s in inlist]
    else:
        return [s.replace(replacefrom, replaceto) for s in inlist]


def getElementFromList(inlist,idx,slc=None):
    '''
    For selecting a particular element or slice from a list 
    within a nipype connect statement.
    If the slice is longer than the list, this returns the list
    '''
    if not slc:
        outlist=inlist[idx]
    else:
        if slc == -1:
            outlist=inlist[idx:]
        else:
            outlist=inlist[idx:slc]
    return outlist


def writeIdentity(**kwargs):
    '''
    Tautology with a dictionary
    Has its uses for meta-interfaces
    '''
    return kwargs# other dictionary functions