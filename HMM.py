"""A generalized Hidden Markov Model that employs the Viterbi algorithm."""
import math
import string
import re
import math
from sys import argv


def DNA(filename):
    """
    Reads in a file and returns a clean string of the file content.
    Parameters
    ----------
    filename : textfile        
    type:String
    Returns
    -------
    out : A string of nucleotide characters {A, C, T, G}
    rtype: String
    """
    with open(filename) as f:
        input_text=f.read()
        x=re.sub('[^acgtACGT]', '',input_text)
    return x.upper()
    infile.close()



def transform(x):
    """
    coverts every intial probability vector to log.
    Parameters
    ----------
    x : initial probability distribution vector of length N         
    
    Returns
    -------
    out : A list of log eqivalent of x
    rtype:float
    """
    for i in range(len(x)):
        x[i]= math.log(x[i])  
    return x


def pmatrix(x):
    """
    computes the values of either the transition or Emission probability matrix to log equivalent
    
    Parameters
    ----------
    x : transition or emission matrix         
    
    Returns
    -------
    out : log matrix
    rtype: list
    """
    for i in range(len(x)):
        for j in range(len(x[0])):
          if x[i][j]==0:
              x[i][j]=-999999999
          else:
              x[i][j]= math.log(x[i][j]) 
    return x


def Transmission(tp):
    """Translates Pmatrix(x) to a dictionary equivalent
    Parameters
    ----------
    tp :  function Pmatrix(x)
    type: list
    
    Returns
    -------
    out : A dictionary of Transition Probability
    rtype: dict
    """
    d={}
    for j in range(len (state)):
        for i in range(len(tp)):
           d[state[i],state[j]] = tp[i][j]
    return d


def Emmission(ep):
    """translates Pmatrix(x) to a dictionary equivalent
    Parameters
    ----------
    ep : function Pmatrix(x)    
    type: list
    Returns
    -------
    out : A dictionary of Emission Probability
    rtype: dict
    """
    d={}
    for j in range(len(e)):
        for i in range(len(state)):
           d[state[i],e[j]] = ep[i][j]
    return d


def Alpha(row,col):
    """ initalise and Alpha Matrix,computes the first colum of the matrix 
    
    Parameters
    ----------
    row : State
    col : function DNA(filename): value return by this function
    
    Returns
    -------
    out : A Matrix with the first column computed as the initailization of Alpha matrix 
    rtype=list
    """
    j=0
    dmatrix=[]
    for i in range(len(state)):
                dmatrix+=[[0]*len(col)]
    df=dmatrix
    for i in range(len(state)):
         df[i][j]=Pi[i]+ EP[(state[i],strng[0])] 
    return df



def maxcolno(df,j,i):
    """ 
    outputs the maximum no and its index from a list of iteration(max(alpha[j]*T[j][i])
    Parameters
    ----------
    df : Alpha matrix
     j :  column
     i : index
    df type: list
    j type: int
    i type:int
    Returns
    -------
    out : x: max no from computation
          t: the index of x
    rtype: x is a float value, t is an int value 
    """
    lst=[]
    for k in range(len(state)):
        lst+=[df[k][j]+TP[(state[k],state[i])]]
    x=max(lst)
    t=lst.index(max(lst))
    return x,t


def maxindex(df,j):
    """ 
    outputs the maximum no and the eqivalent index of the column
    Parameters
    ----------
    df : Alpha matrix
     j :  column
    type: df is a list, j is an int
    Returns
    -------
    out : x: max no 
          w: the index of x
    type: df is a list, j is an int
    """
    maxn=[]
    for i in range(len(state)): 
          maxn+=[df[i][j]]
    x=max(maxn) 
    w= maxn.index(x)
    return x,w
    


def iteration(df):
    """ 
    iterates through the length of the Alpha matrix excluding the first column
    Parameters
    ----------
    df : Alpha matrix
    Returns
    -------
    out : df: A complete computed matrix of HMM
          vert; A matrix holding the indicies of max computed from a list of values 
    rtype: list
    """
    vert=[]
    for i in range(len(state)):
        vert+=[[0]* len(strng)]
               
    for j in range(1,len(strng)):
        for i in range(len(state)):
            df[i][j]=EP[(state[i],strng[j])]+maxcolno(df,j-1,i)[0]
            vert[i][j]=maxcolno(df,j-1,i)[1]
                     
    return df,vert



def pathv(df,vert):
    """ 
    Do a Backtrace of hidden path
    Parameters
    ----------
    df : Alpha matrix
    vert: Matrix holding indicies of states
    Returns
    -------
    out : string of best possible path
    rtype: string
    """
    finalpath=''
    j=len(strng)-1
    maxlastcol= maxindex(df,j)[0]
    value=maxindex(df,j)[1]
    bpath=state[value]
    i=value
    for j in range(len(strng)-1,0,-1):
        value=vert[i][j]                
        statevalue=state[value]
        finalpath=statevalue+","+finalpath
        i=value    
    return finalpath+bpath,"The best path probability=", maxlastcol
    


if __name__ == '__main__':

    state=['S1','S2','S3']
    pi=[.3,.2,.5]
    Pi=transform(pi)
    e=['A','C','T','G']
    tpmatrix=[[.5,0,.5],[.25,.5,.25],[.1,.4,.5]]
    tp=pmatrix(tpmatrix)
    TP= Transmission(tp)
    epmatrix=[[.3,.1,.4,.2],[.1,.5,.1,.3],[.25,.25,.25,.25]]
    ep=pmatrix(epmatrix)
    EP=Emmission(epmatrix)
    #myfile=argv[1]
    myfile="seq.fasta.txt"
    strng=DNA(myfile)
    df=Alpha(state,strng)
    HMM=iteration(df)
    x=HMM[0]
    y=HMM[1]
    print pathv(x,y)
    

