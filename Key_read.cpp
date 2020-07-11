#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>

using namespace std;

static vector<int> private_key;
static vector<vector<int> > U;
static vector<vector<int> > HybridMatrix;
static int nblock_bio;
static int nblock_pri;

static vector<int> HYBRID_CODE;

void multiplication(vector<vector<int> > permMatrix)
{
    HybridMatrix.resize(U.size()); //declaring no. of rows to be same
    int i, j, k;
    for (i = 0; i < U.size(); i++)
    {
        for (j = 0; j < U.size(); j++)
        {
            HybridMatrix[i].push_back(0);
            for (k = 0; k < U.size(); k++)
                HybridMatrix[i][j] += U[i][k]*permMatrix[k][j];
        }
    }
}



void permutation_and_multiplication()
{
    //6 permutation each 3x3 size
    int p0[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    int p1[3][3]={{1,0,0},{0,0,1},{0,1,0}};
    int p2[3][3]={{0,1,0},{1,0,0},{0,0,1}};
    int p3[3][3]={{0,0,1},{1,0,0},{0,1,0}};
    int p4[3][3]={{0,1,0},{0,0,1},{1,0,0}};
    int p5[3][3]={{0,0,1},{0,1,0},{1,0,0}};

    int row_covered=0;
    vector<vector<int> > permMatrix( U.size() , vector<int> (U.size(), 0));  // declare a matrix of same dimension as U, and initialize all elements as 0

    for(int i=0;i<private_key.size();i++)
    {
        int perm_select=private_key[i]%6; //to select one of the above three permutations
        if(row_covered<U.size())
        {
            int start=3*i;
            int endd=(3*i)+3;
            for(int j=start;j<endd;j++)
            {
                for(int k=start;k<endd;k++)
                {
                    if(perm_select==0)
                        permMatrix[j][k]=p0[j-(3*i)][k-(3*i)];
                    else if(perm_select==1)
                        permMatrix[j][k]=p1[j-(3*i)][k-(3*i)];
                    else if(perm_select==2)
                        permMatrix[j][k]=p2[j-(3*i)][k-(3*i)];
                    else if(perm_select==3)
                        permMatrix[j][k]=p3[j-(3*i)][k-(3*i)];
                    else if(perm_select==4)
                        permMatrix[j][k]=p5[j-(3*i)][k-(3*i)];
                    else if(perm_select==5)
                        permMatrix[j][k]=p5[j-(3*i)][k-(3*i)];
                }
            }
            row_covered+=3; // 3 rows covered
        }
        else
            break;
    }

    multiplication(permMatrix);
}



void addExtra(int Diff)
{
    if(Diff<0)//extra row
    {
        for(int j=0;j<abs(Diff);j++)
        {
            int col_count=0;
            int col_size=U[0].size(); //size of col //assuming same no. of elements

            U.resize(U.size()+1); //adding 1 more row
            int last_row=U.size()-1; // last row index
            for(int i=0;i<private_key.size();i++)
            {
                if(i%2==0)
                {
                    if(col_count<col_size)
                    {
                        int index=private_key[i];
                        int val=U[j][index];
                        U[last_row].push_back(val);
                        col_count++;
                    }
                }
                else
                {
                    if(col_count<col_size)
                    {
                        int index=private_key[i];
                        int val=U[index][j];
                        U[last_row].push_back(val);
                        col_count++;
                    }
                }
            }
        }

    }
    else if(Diff>0) //extra col
    {
        for(int j=0;j<abs(Diff);j++)
        {
            int row_count=0;
            int col_size=U[0].size(); // col size of matrix U
            for(int i=0;i<private_key.size();i++)
            {
                if(i%2==0) // if even, pick value from first row
                {
                    if(row_count<=U.size())
                    {
                        int index=private_key[i];
                        int val=U[j][index];
                        U[i].push_back(val);
                        row_count++;
                    }
                }
                else
                {
                    if(row_count<=U.size())
                    {
                        int index=private_key[i];
                        int val=U[index][j];
                        U[i].push_back(val);
                        row_count++;
                    }
                }
            }
        }

    }

}



vector<vector<int> > matrix_generation(int first_size, int second_size,vector<int> first, vector<int> second, int block1,int block2)
{
    vector<vector<int> > Matrix(block1+block2); // max rows that can be formed right now : block1 + block2
    int *pointFirst=first.data(); //points to the first element of first vector
    int *pointSecond=second.data(); //points to the second element of first vector
    int nxt_row_to_fill=0; // for the matrix


    int actual_block_1 = block1,actual_block_2 = block2; // actual values of block

    int frst_cmp_size=first.size(); // for compairing the no. of items/elements/components
    int snd_cmp_size=second.size();

    for(int i=0;i<private_key.size();i++) //run full private key vector
    {
        if(block1!=0 && block2!=0)
        {
            if(i%2==0) // even index, work with first vector
            {
                int item = ceil(first_size/actual_block_1); // no. of elements we need to add to each row when using this vector
                int blc_to_use = private_key[i]; // no. of blocks that we have to add continuously

                if(blc_to_use<=block1)
                {
                    int row_range=nxt_row_to_fill+blc_to_use;
                    for(int j=nxt_row_to_fill;j<row_range;j++)
                    {
                        if(item<=frst_cmp_size)
                        {
                            for(int k=0;k<item;k++)
                            {
                                Matrix[j].push_back(*pointFirst++);
                            }
                            frst_cmp_size=frst_cmp_size-item;
                        }
                        else if(item>frst_cmp_size)
                        {
                            for(int k=0;k<frst_cmp_size;k++)
                            {
                                Matrix[j].push_back(*pointFirst++);
                            }
                            frst_cmp_size=0;
                        }
                    }
                    nxt_row_to_fill=nxt_row_to_fill+blc_to_use;
                    block1=block1-blc_to_use; // no. of blocks left to use
                }
                else if(blc_to_use>block1)
                {
                    int row_range=nxt_row_to_fill+block1;
                    for(int j=nxt_row_to_fill;j<row_range;j++)
                    {
                        if(item<=frst_cmp_size)
                        {
                            for(int k=0;k<item;k++)
                            {
                                Matrix[j].push_back(*pointFirst++);
                            }
                            frst_cmp_size=frst_cmp_size-item;
                        }
                        else if(item>frst_cmp_size)
                        {
                            for(int k=0;k<frst_cmp_size;k++)
                            {
                                Matrix[j].push_back(*pointFirst++);
                            }
                            frst_cmp_size=0;
                        }
                    }
                    nxt_row_to_fill=nxt_row_to_fill+blc_to_use;
                    block1=0; // no. of blocks left to use
                }
            }
            else if(i%2==1)
            {
                int item = ceil(second_size/actual_block_2);
                int blc_to_use = private_key[i];

                if(blc_to_use<=block2)
                {
                    int row_range=nxt_row_to_fill+blc_to_use;
                    for(int j=nxt_row_to_fill;j<row_range;j++)
                    {
                        if(item<=snd_cmp_size)
                        {
                            for(int k=0;k<item;k++)
                            {
                                Matrix[j].push_back(*pointSecond++);
                            }
                            snd_cmp_size=snd_cmp_size-item;
                        }
                        else if(item>snd_cmp_size)
                        {
                            for(int k=0;k<snd_cmp_size;k++)
                            {
                                Matrix[j].push_back(*pointSecond++);
                            }
                            snd_cmp_size=0;
                        }
                    }
                    nxt_row_to_fill=nxt_row_to_fill+blc_to_use;
                    block2=block2-blc_to_use; // no. of blocks left to use
                }
                else if(blc_to_use>block2)
                {
                    int row_range=nxt_row_to_fill+block2;
                    for(int j=nxt_row_to_fill;j<row_range;j++)
                    {
                        if(item<=snd_cmp_size)
                        {
                            for(int k=0;k<item;k++)
                            {
                                Matrix[j].push_back(*pointSecond++);
                            }
                            snd_cmp_size=snd_cmp_size-item;
                        }
                        else if(item>snd_cmp_size)
                        {
                            for(int k=0;k<snd_cmp_size;k++)
                            {
                                Matrix[j].push_back(*pointSecond++);
                            }
                            snd_cmp_size=0;
                        }
                    }
                    nxt_row_to_fill=nxt_row_to_fill+block2;
                    block2=0;
                }
            }
        }
        else if(block1!=0)
        {
            int item = ceil(first_size/actual_block_1);
            int row_range=nxt_row_to_fill+block1;
            for(int j=nxt_row_to_fill;j<row_range;j++)
            {
                if(item<=frst_cmp_size)
                {
                    for(int k=0;k<item;k++)
                    {
                        Matrix[j].push_back(*pointFirst++);
                    }
                    frst_cmp_size=frst_cmp_size-item;
                }
                else if(item>frst_cmp_size)
                {
                    for(int k=0;k<frst_cmp_size;k++)
                    {
                        Matrix[j].push_back(*pointFirst++);
                    }
                    frst_cmp_size=0;
                }
            }
            nxt_row_to_fill=nxt_row_to_fill+block1;
            block1=0; // no. of blocks left to use
            break;
        }
        else if(block2!=0)
        {
            int item = ceil(second_size/actual_block_2);

            int row_range=nxt_row_to_fill+block2;
            for(int j=nxt_row_to_fill;j<row_range;j++)
            {
                if(item<=snd_cmp_size)
                {
                    for(int k=0;k<item;k++)
                    {
                        Matrix[j].push_back(*pointSecond++);
                    }
                    snd_cmp_size=snd_cmp_size-item;
                }
                else if(item>snd_cmp_size)
                {
                    for(int k=0;k<snd_cmp_size;k++)
                    {
                        Matrix[j].push_back(*pointSecond++);
                    }
                    snd_cmp_size=0;
                }
            }
            nxt_row_to_fill=nxt_row_to_fill+block2;
            block2=0; // no. of blocks left to use
            break;
        }
    }
    return Matrix;
}


vector<int> padd_fill(vector<int> key,int index,int old_size)
{
    for(int i=old_size;i<key.size();i++)
    {
        key[i]=key[index];
        index=index+private_key[i-old_size+1];
    }

    return key;
}

vector<int> key_to_vector_generator(char* file_name)
{
    vector<int> key;//biometric vector
    fstream face_key(file_name,ios::in); //opening biometric key consisting file n read mode
    if(!face_key)
        exit(1);
    else
    {
        while(!face_key.eof())
        {
            string s;
            face_key>>s;

            if(s!=" " || s!="\n")
            {
                stringstream val(s);
                int value=0;
                val>>value;
                key.push_back(value); // creating a vector with the bio keys
            }
        }
        face_key.close();
    }
    return key;
}


int main()
{
    char* bio_name="face_key.txt";
    char* pri_name="prime_key.txt";
    char* private_name="private_key.txt";

    vector<int> bio_key=key_to_vector_generator(bio_name);//biometric vector
    vector<int> pri_key=key_to_vector_generator(pri_name);//prime key vector
    private_key=key_to_vector_generator(private_name);//prime key vector

    long int m=bio_key.size();
    long int n=pri_key.size();
    long int s = m +  n; // s = m + n contains both their size
    long int q = ceil(sqrt(s));

    int nz1=q-(m%q); // padding for bio_key vector
    int nz2=q-(n%q); // padding for pri_key vector

    bio_key.resize(m+nz1,-1); //resize bio vector to accomodate padding, initially -1, as pixels and key can't have -ve values
    pri_key.resize(n+nz2,-1); //resize prime vector to accomodate padding

    int index=private_key[0];
    bio_key=padd_fill(bio_key,index,m); //padding fill
    pri_key=padd_fill(pri_key,index,n);


    nblock_bio = bio_key.size()/q; //dividing vector into blocks
    nblock_pri = pri_key.size()/q;

    int Pad=nz1+nz2;
    int old_bio_size=bio_key.size();
    int old_pri_size=pri_key.size();

    bio_key.resize(bio_key.size()+Pad,-1); // adding the padded size to the two vectors
    pri_key.resize(pri_key.size()+Pad,-1);

    bio_key=padd_fill(bio_key,index,old_bio_size);
    pri_key=padd_fill(pri_key,index,old_pri_size);

    // now which vector will start the matrix 'U' is decided by first component of private key module 2
    index = private_key[0]%2;
    if(index==0)
      U= matrix_generation(bio_key.size(),pri_key.size(),bio_key,pri_key,nblock_bio,nblock_pri);
    else
      U= matrix_generation(pri_key.size(),bio_key.size(),pri_key,bio_key,nblock_pri,nblock_bio);


    //little flaw in the maths...which is calculated, hence change of algo
    //int padTotal=(q*q)-(m+n);
    //int Diff=Pad-padTotal;
    int Diff=U.size()-U[0].size(); //assuming that all line have same no. of columns (row-col)
    addExtra(Diff);

    // this here makes the matrix dimension divisible by 3
    int row=U.size();
    int col=U[0].size();
    if(row%3!=0)
    {
        int time=row%3;
        addExtra(-1*time);
    }
    if(col%3!=0)
    {
        int time=col%3;
        addExtra(time);
    }

    permutation_and_multiplication(); //for permutation and post-multiplication of U and permutedMatrix

    for(int i=0;i<HybridMatrix.size();i++)
    {
        for(int j=0;j<HybridMatrix[i].size();j++)
        {
            HYBRID_CODE.push_back(HybridMatrix[i][j]);
        }
    }

    ofstream hyb("hybrid_code.txt");
    for(int i=0;i<HYBRID_CODE.size();i++)
    {
        hyb<<HYBRID_CODE[i];
    }
    hyb.close();

    for(int i=0;i<HYBRID_CODE.size();i++)
    {
        cout<<HYBRID_CODE[i];
    }

    return 0;
}

