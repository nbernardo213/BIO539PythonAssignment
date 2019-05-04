import pandas as pd
import plotnine as p9
import sys

def read_dat(Sequence,k=0):
    """
    Loads a sequence and sets k-mer length to k.
    k must be less than total sequence length and must be an integer.
    Sequence should only contain letters (ATCG).

    Sequence(String): String of characters representing a sequence
    k(int): Length of desired k-mers

    Returns:
        Unique_len(int): Number of observed k-mers
    """
    k=int(k)
    test_len=len(Sequence)
    for i in range(test_len - k + 1):
        new=Sequence[i:k]
        if i==0:
            t=k
            unique=[]
            unique.append(new)
        else:
            p=0
            for j in range(len(unique)):
                if unique[j] == new:
                    p=1
                else:
                    p=p
            if p==0:
                unique.append(new)
        k=k+1
    Unique_len=len(unique)
    print("The number of observed k-mers is:")
    print(Unique_len)
    return(Unique_len)

def gen_poss(Sequence_Length,t=0):
    """
    Calculates the possible number of k-mers given a sequence.

    Sequence_Length (int): length of the original sequence.
    t (int): desirerd k-mer length

    returns:
        possible (int): the total possible k-mers in the sequence of length t
    """
    t=int(t)
    if Sequence_Length-t+1>4**t:
        possible = 4**t
    else:
        possible = Sequence_Length-t+1
    print("The number of possible k-mers is:")
    print(possible)
    return(possible)

def kmer_df(Sequence):
    """
    Creates a dataframe from a sequence containing all possible lengths k and
    all possible and observed k-mer counts.

    Sequence (str): sequence in the format of a string (i.e., ATTCGA)

    Returns:
        df (dataframe): dataframe containing observed and possible number of
                        k-mers for all possible k-mer lengths
    """
    poss_k=[]
    poss_possible=[]
    poss_unique=[]
    test_len=len(Sequence)
    test=Sequence
    for m in range(test_len+1):
        poss_k.append(m)
        for i in range(test_len-m+1):
            new=test[i:m]
            if i==0:
                t=m
                unique=[]
                unique.append(new)
            else:
                p=0
                for j in range(len(unique)):
                    if unique[j] == new:
                        p=1
                    else:
                        p=p
                if p==0:
                    unique.append(new)
            m=m+1

            if i==0:
                possible = 0
                if test_len-t+1>4**t:
                    possible = 4**t
                else:
                    possible = test_len-t+1
        m=t
        poss_possible.append(possible)
        poss_unique.append(len(unique))

    data=[poss_unique,poss_possible]
    df=pd.DataFrame(data)
    df.loc[2]=df.loc[0]/df.loc[1]
    df.loc[3]=poss_k
    df = df.rename(index={0:'Observed',1:'Possible',2:'Complexity',3:'k'})
    df=df.T
    df
    print(df)
    return(df)

def complexity(df):
    """
    generates the complexity of a sequence.

    df (dataframe): dataframe containing observed and possible k-mers
                    for all possible lengths

    returns:
        comp (double): calculated complexity of the sequence (observed/possible)
    """
    comp=df['Observed'].sum()/df['Possible'].sum()
    print("The complexity of the sequence is:")
    print(comp)

def plot_df(df):
    """
    Generates a plot of complexity of sequence over k-mer length k.

    df (dataframe): dataframe containing observed and possible k-mers for a
                    sequence of characters.

    returns:
        plot (plot): graph of complexity over length k.
    """
    return(p9.ggplot(data=df, mapping=p9.aes(x='k',y='Complexity'))+p9.geom_point()+p9.geom_line())

if __name__ == '__main__':
    if len(sys.argv)==3:
        myfile = sys.argv[1]
        lengthk = sys.argv[2]
        lengthk=int(lengthk)
        u=0
        for z in range(len(myfile)):
            if myfile[z]!= 'A' and myfile[z]!='T' and myfile[z]!='C' and myfile[z]!='G':
                u=1
            else:
                u=u
        if u==1:
            print("Sequence can only contain the letters: A T C G")
        else:
            if lengthk>len(myfile):
                print("Length of k must be smaller than the length of the sequence!")
            else:
                seq=read_dat(myfile,lengthk)
                possible=gen_poss(seq,lengthk)

            all_df=kmer_df(myfile)
            complex=complexity(all_df)
            plot=plot_df(all_df)
            print(plot)
    else:
        if len(sys.argv)>3:
            print("Too many arguments on command line.")
        else:
            if len(sys.argv)==2:
                myfile = sys.argv[1]
                all_df=kmer_df(myfile)
                complex=complexity(all_df)
                plot=plot_df(all_df)
                print(plot)
            else:
                print("Please enter a sequence, such as ATTGCA, followed by a space and then a desired length of k-mer (k-mer length optional). An example of inputs appears as: ATTCGA 3")
