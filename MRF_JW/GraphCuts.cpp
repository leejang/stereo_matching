/*
 * File: GraphCuts.cpp
 *
 * This file originzted from the MRF 2.2 Library
 * slightly modfied by Jangwon Lee, 11/2016 
 * Email: leejang@indiana.edu
 */

#include "energy.h"
#include "graph.h"
#include "GraphCuts.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "string.h"
#define MAX_INTT 1000000000



/**************************************************************************************/

void GraphCuts::initialize_memory()
{

    m_lookupPixVar = (PixelType *) new PixelType[m_nPixels];
    m_labelTable   = (LabelType *) new LabelType[m_nLabels];

    terminateOnError( !m_lookupPixVar || !m_labelTable ,"Not enough memory");

    for ( int i = 0; i < m_nLabels; i++ )
        m_labelTable[i] = i;
}

/**************************************************************************************/

void GraphCuts::commonGridInitialization(PixelType width, PixelType height, int nLabels)
{
    terminateOnError( (width < 0) || (height <0) || (nLabels <0 ),"Illegal negative parameters");

    m_width       = width;
    m_height      = height;
    m_nPixels     = width*height;
    m_nLabels     = nLabels;
    m_grid_graph  = 1;
        

}

/**************************************************************************************/

void GraphCuts::setParameters(int numParam, void *param)
{
    if (numParam != 1 ) 
        printf("\nInvalid number of parameters, can only set one parameter\\that is boolean label order\n");
    else
    {
        m_random_label_order = *((bool *) param);
    }
}

/**************************************************************************************/

void GraphCuts::commonNonGridInitialization(PixelType nupixels, int num_labels)
{
    terminateOnError( (nupixels <0) || (num_labels <0 ),"Illegal negative parameters");

    m_nLabels        = num_labels;
    m_nPixels        = nupixels;
    m_grid_graph         = 0;

    m_neighbors = (LinkedBlockList *) new LinkedBlockList[nupixels];

    terminateOnError(!m_neighbors,"Not enough memory");

}

/**************************************************************************************/

void GraphCuts::commonInitialization()
{
    m_needToFreeV        = 0;
    m_random_label_order = 1;
    initialize_memory();
}

/**************************************************************************************/
/* Use this constructor only for grid graphs                                          */
GraphCuts::GraphCuts(PixelType width,PixelType height,LabelType nLabels,EnergyFunction *eng):MRF(width,height,nLabels,eng)
{
    commonGridInitialization(width,height,nLabels);
    
    m_labeling           = (LabelType *) new LabelType[m_nPixels];
    terminateOnError(!m_labeling,"out of memory");
    for ( int i = 0; i < m_nPixels; i++ ) m_labeling[i] = (LabelType) 0;

    commonInitialization(); 
}

/**************************************************************************************/
/* Use this constructor for general graphs                                            */
GraphCuts::GraphCuts(PixelType nPixels,int nLabels,EnergyFunction *eng):MRF(nPixels,nLabels,eng) 
{

    commonNonGridInitialization(nPixels, nLabels);
    
    m_labeling           = (LabelType *) new LabelType[m_nPixels];
    terminateOnError(!m_labeling,"out of memory");
    for ( int i = 0; i < nPixels; i++ ) m_labeling[i] = (LabelType) 0;

    commonInitialization(); 
}

/**************************************************************************************/

void GraphCuts::setData(EnergyTermType* dataArray)
{
    m_datacost  = dataArray;
}


/**************************************************************************************/

// Truncated L1
void GraphCuts::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
{
    int i, j;
    CostVal cost;

    m_needToFreeV = 1;

    m_smoothcost = (CostVal *) new CostVal[m_nLabels*m_nLabels];
    if (!m_smoothcost) { fprintf(stderr, "Not enough memory!\n"); exit(1); }


    for (i=0; i<m_nLabels; i++)
        for (j=i; j<m_nLabels; j++)
        {
            cost = (CostVal)(j - i);
            if (cost > smoothMax) cost = smoothMax;
            m_smoothcost[i*m_nLabels + j] = m_smoothcost[j*m_nLabels + i] = cost*lambda;
        }

}

#if 0
// Pott Model
void GraphCuts::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
{
    int i, j;
    CostVal cost;

    m_needToFreeV = 1;

    m_smoothcost = (CostVal *) new CostVal[m_nLabels*m_nLabels];
    if (!m_smoothcost) { fprintf(stderr, "Not enough memory!\n"); exit(1); }


    for (i=0; i<m_nLabels; i++)
        for (j=i; j<m_nLabels; j++) {

            if (j == i)
              cost = (CostVal) 0;
            else
              cost = (CostVal) 10;

            m_smoothcost[i*m_nLabels + j] = m_smoothcost[j*m_nLabels + i] = cost;
        }

}
#endif

/**************************************************************************************/

void GraphCuts::setSmoothness(SmoothCostGeneralFn cost)
{
    m_smoothFnPix   = cost;
}

/**************************************************************************************/

GraphCuts::EnergyType GraphCuts::dataEnergy()
{
    
    if (m_dataType == ARRAY) 
        return(giveDataEnergyArray());
    else
        fprintf(stderr, "Did not support Function Type Data Energy!\n");

    return(0);
}

/**************************************************************************************/
    
GraphCuts::EnergyType GraphCuts::giveDataEnergyArray()
{
    EnergyType eng = (EnergyType) 0;


    for ( int i = 0; i < m_nPixels; i++ )
        eng = eng + m_datacost(i,m_labeling[i]);

    return(eng);
}


/**************************************************************************************/

GraphCuts::EnergyType GraphCuts::smoothnessEnergy()
{
    
    if (m_grid_graph)
    {
        if ( m_smoothType != FUNCTION )
        {
            if (m_varWeights) return(fprintf(stderr, "Did not support weighted array for data cost!\n"));
            else return(giveSmoothEnergy_G_ARRAY());
        }
        else return(fprintf(stderr, "Did not support Function Type Data Energy!\n"));
    }
    else
    {
        fprintf(stderr, "Did not support Non Grid Type!\n");
    }
    return(0);

}

/**************************************************************************************/

GraphCuts::EnergyType GraphCuts::giveSmoothEnergy_G_ARRAY()
{

    EnergyType eng = (EnergyType) 0;
    int x,y,pix;


    for ( y = 0; y < m_height; y++ )
        for ( x = 1; x < m_width; x++ )
        {
            pix = x+y*m_width;
            eng = eng + m_smoothcost(m_labeling[pix],m_labeling[pix-1]);
        }

    for ( y = 1; y < m_height; y++ )
        for ( x = 0; x < m_width; x++ )
        {
            pix = x+y*m_width;
            eng = eng + m_smoothcost(m_labeling[pix],m_labeling[pix-m_width]);
        }

    return(eng);
    
}

/**************************************************************************************/

void GraphCuts::setLabelOrder(bool RANDOM_LABEL_ORDER)
{
    m_random_label_order = RANDOM_LABEL_ORDER;
}

/****************************************************************************/
/* This procedure checks if an error has occured, terminates program if yes */

void GraphCuts::terminateOnError(bool error_condition,const char *message)

{ 
   if  (error_condition) 
   {
      fprintf(stderr, "\n %s \n", message);
      exit(1);
    }
}

/****************************************************************************/

void GraphCuts::clearAnswer()
{
    if (!m_labeling ) {printf("First initialize algorithm" );exit(0);}
    memset(m_labeling, 0, m_nPixels*sizeof(Label));
}

/****************************************************************************/

void GraphCuts::scramble_label_table()
{
   LabelType r1,r2,temp;
   int num_times,cnt;

   num_times = m_nLabels*2;
   srand(clock());

   for ( cnt = 0; cnt < num_times; cnt++ )
   {
      r1 = rand()%m_nLabels;  
      r2 = rand()%m_nLabels;  

      temp             = m_labelTable[r1];
      m_labelTable[r1] = m_labelTable[r2];
      m_labelTable[r2] = temp;
   }
}

/**************************************************************************************/

void GraphCuts::add_t_links_ARRAY(Energy *e,Energy::Var *variables,int size,LabelType alpha_label)
{
    for ( int i = 0; i < size; i++ )
        e -> add_term1(variables[i], m_datacost(m_lookupPixVar[i],alpha_label),
                                     m_datacost(m_lookupPixVar[i],m_labeling[m_lookupPixVar[i]]));

}

/**************************************************************************************/

void GraphCuts::add_t_links_FnPix(Energy *e,Energy::Var *variables,int size,LabelType alpha_label)
{
    for ( int i = 0; i < size; i++ )
        e -> add_term1(variables[i], m_dataFnPix(m_lookupPixVar[i],alpha_label),
                                     m_dataFnPix(m_lookupPixVar[i],m_labeling[m_lookupPixVar[i]]));

}
/**************************************************************************************/

void GraphCuts::setNeighbors(PixelType pixel1, int pixel2, EnergyTermType weight)
{

    assert(pixel1 < m_nPixels && pixel1 >= 0 && pixel2 < m_nPixels && pixel2 >= 0);
    assert(m_grid_graph == 0);

    Neighbor *temp1 = (Neighbor *) new Neighbor;
    Neighbor *temp2 = (Neighbor *) new Neighbor;

    temp1->weight  = weight;
    temp1->to_node = pixel2;

    temp2->weight  = weight;
    temp2->to_node = pixel1;

    m_neighbors[pixel1].addFront(temp1);
    m_neighbors[pixel2].addFront(temp2);
    
}

/**************************************************************************************/

GraphCuts::~GraphCuts()
{


    delete [] m_labeling;
    if ( ! m_grid_graph ) delete [] m_neighbors;            
    delete [] m_labelTable;
    delete [] m_lookupPixVar;
    if (m_needToFreeV) delete [] m_smoothcost;
}

/**************************************************************************************/

Swap::Swap(PixelType width,PixelType height,int num_labels, EnergyFunction *eng):GraphCuts(width,height,num_labels,eng)
{
    m_pixels = new PixelType[m_nPixels];
    terminateOnError( !m_pixels ,"Not enough memory");
}

/**************************************************************************************/

Swap::~Swap()
{
    delete [] m_pixels;
}

/**************************************************************************************/

GraphCuts::EnergyType Swap::swap(int max_num_iterations)
{
    return(start_swap(max_num_iterations)); 
}

/**************************************************************************************/


GraphCuts::EnergyType Swap::start_swap(int max_num_iterations )
{
    
    int curr_cycle = 1;
    EnergyType new_energy,old_energy;
    

    new_energy = dataEnergy()+smoothnessEnergy();

    //old_energy = new_energy+1; // this doesn't work for large float energies
    old_energy = -1;

    while (old_energy < 0 || (old_energy > new_energy  && curr_cycle <= max_num_iterations))
    {
        old_energy = new_energy;
        new_energy = oneSwapIteration();
        
        curr_cycle++;   
    }
    //printf(" swap energy %d",new_energy);
    return(new_energy);
}

/**************************************************************************************/

GraphCuts::EnergyType Swap::oneSwapIteration()
{
    int next,next1;
   
    if (m_random_label_order) scramble_label_table();
        

    for (next = 0;  next < m_nLabels;  next++ )
        for (next1 = m_nLabels - 1;  next1 >= 0;  next1-- )
            if ( m_labelTable[next] < m_labelTable[next1] )
            {
                perform_alpha_beta_swap(m_labelTable[next],m_labelTable[next1]);
            }

    return(dataEnergy()+smoothnessEnergy());
}

/**************************************************************************************/

GraphCuts::EnergyType Swap::alpha_beta_swap(LabelType alpha_label, LabelType beta_label)
{
    terminateOnError( alpha_label < 0 || alpha_label >= m_nLabels || beta_label < 0 || beta_label >= m_nLabels,
        "Illegal Label to Expand On");
    perform_alpha_beta_swap(alpha_label,beta_label);
    return(dataEnergy()+smoothnessEnergy());
}
/**************************************************************************************/

void Swap::add_t_links_ARRAY_swap(Energy *e,Energy::Var *variables,int size,
                                            LabelType alpha_label, LabelType beta_label,
                                            PixelType *pixels)
{
    for ( int i = 0; i < size; i++ )
        e -> add_term1(variables[i], m_datacost(pixels[i],alpha_label),
                                     m_datacost(pixels[i],beta_label));

}
    
/**************************************************************************************/

void Swap::add_t_links_FnPix_swap(Energy *e,Energy::Var *variables,int size,
                                            LabelType alpha_label, LabelType beta_label,
                                            PixelType *pixels)
{
    for ( int i = 0; i < size; i++ )
        e -> add_term1(variables[i], m_dataFnPix(pixels[i],alpha_label),
                                     m_dataFnPix(pixels[i],beta_label));

}



/**************************************************************************************/

void Swap::perform_alpha_beta_swap(LabelType alpha_label, LabelType beta_label)
{
    PixelType i,size = 0;
    Energy *e = new Energy();


    for ( i = 0; i < m_nPixels; i++ )
    {
        if ( m_labeling[i] == alpha_label || m_labeling[i] == beta_label)
        {
            m_pixels[size]    = i;
            m_lookupPixVar[i] = size;
            size++;
        }
    }

    if ( size == 0 ) return;


    Energy::Var *variables = (Energy::Var *) new Energy::Var[size];
    if (!variables) { fprintf(stderr, "Not enough memory!\n"); exit(1); }
    

    for ( i = 0; i < size; i++ )
        variables[i] = e ->add_variable();

    if (m_dataType == ARRAY)
        add_t_links_ARRAY_swap(e,variables,size,alpha_label,beta_label,m_pixels);
    else 
        fprintf(stderr, "Did not support Function Type Data Energy!\n");

    if  (m_grid_graph) {
        if ( m_smoothType != FUNCTION ) {
            if (m_varWeights)
                fprintf(stderr, "Did not support weighted array for data cost!\n");
            else
                set_up_swap_energy_G_ARRAY(size,alpha_label,beta_label,m_pixels,e,variables);
        } else {
            fprintf(stderr, "Did not support Function Type Data Energy!\n");
        }
        
    }
    else
    {
        fprintf(stderr, "Did not support Non Grid Type!\n");
    }
        

    e -> minimize();

    for ( i = 0; i < size; i++ )
        if ( e->get_var(variables[i]) == 0 )
            m_labeling[m_pixels[i]] = alpha_label;
        else m_labeling[m_pixels[i]] = beta_label;


    delete [] variables;
    delete e;

}

/**************************************************************************************/

void Swap::set_up_swap_energy_G_ARRAY(int size,LabelType alpha_label,LabelType beta_label,
                                               PixelType *pixels,Energy* e, Energy::Var *variables)

{
    PixelType nPix,pix,i,x,y;
    


    for ( i = 0; i < size; i++ )
    {
        pix = pixels[i];
        y = pix/m_width;
        x = pix - y*m_width;

        if ( x > 0 )
        {
            nPix = pix - 1;
    
    
            if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label),
                              m_smoothcost(alpha_label,beta_label),
                              m_smoothcost(beta_label,alpha_label),
                              m_smoothcost(beta_label,beta_label) );
    
                else
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
                                           m_smoothcost(beta_label,m_labeling[nPix]));
    
        }   
        if ( y > 0 )
        {
            nPix = pix - m_width;

            if ( m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label)
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label),
                              m_smoothcost(alpha_label,beta_label),
                              m_smoothcost(beta_label,alpha_label),
                              m_smoothcost(beta_label,beta_label) );
    
                else
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
                                           m_smoothcost(beta_label,m_labeling[nPix]));
        }   

        if ( x < m_width - 1 )
        {
            nPix = pix + 1;
    
            if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label) )
                    e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
                                               m_smoothcost(beta_label,m_labeling[nPix]));
        }   

        if ( y < m_height - 1 )
        {
            nPix = pix + m_width;
    
            if ( !(m_labeling[nPix] == alpha_label || m_labeling[nPix] == beta_label))
                e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
                                           m_smoothcost(beta_label,m_labeling[nPix]));

        }
    }
}

/**************************************************************************************/

void Swap::optimizeAlg(int nIterations)
{
    swap(nIterations);
}

/**************************************************************************************/


GraphCuts::EnergyType Expansion::expansion(int max_num_iterations)
{
    return(start_expansion(max_num_iterations)); 
}

/**************************************************************************************/


GraphCuts::EnergyType Expansion::start_expansion(int max_num_iterations )
{
    
    int curr_cycle = 1;
    EnergyType new_energy,old_energy;
    

    new_energy = dataEnergy()+smoothnessEnergy();


    //old_energy = new_energy+1; // this doesn't work for large float energies
    old_energy = -1;

    while (old_energy < 0 || (old_energy > new_energy  && curr_cycle <= max_num_iterations))
    {
        old_energy = new_energy;
        new_energy = oneExpansionIteration();
        
        curr_cycle++;   
    }

    //printf(" Exp energy %d",new_energy);
    return(new_energy);
}


/**************************************************************************************/

GraphCuts::EnergyType Expansion::alpha_expansion(LabelType label)
{
    terminateOnError( label < 0 || label >= m_nLabels,"Illegal Label to Expand On");

    perform_alpha_expansion(label);
    return(dataEnergy()+smoothnessEnergy());
}


/**************************************************************************************/

void Expansion::perform_alpha_expansion(LabelType alpha_label)
{
    PixelType i,size = 0; 
    Energy *e = new Energy();
    

    
    for ( i = 0; i < m_nPixels; i++ )
    {
        if ( m_labeling[i] != alpha_label )
        {
            m_lookupPixVar[size] = i;
            size++;
        }
    }

        
    if ( size > 0 ) 
    {
        Energy::Var *variables = (Energy::Var *) new Energy::Var[size];
        if ( !variables) {printf("\nOut of memory, exiting");exit(1);}

        for ( i = 0; i < size; i++ )
            variables[i] = e ->add_variable();

        if ( m_dataType == ARRAY ) add_t_links_ARRAY(e,variables,size,alpha_label);
        else  add_t_links_FnPix(e,variables,size,alpha_label);


        if (m_grid_graph) {
            if ( m_smoothType != FUNCTION ) {
                if (m_varWeights)
                    fprintf(stderr, "Did not support weighted array for data cost!\n");
                else
                    set_up_expansion_energy_G_ARRAY(size,alpha_label,e,variables);
            } else {
               fprintf(stderr, "Did not support Function Type Data Energy!\n");
            }
        } else {
            fprintf(stderr, "Did not support Non Grid Type!\n");
        }
        
        e -> minimize();
    
        for ( i = 0,size = 0; i < m_nPixels; i++ )
        {
            if ( m_labeling[i] != alpha_label )
            {
                if ( e->get_var(variables[size]) == 0 )
                    m_labeling[i] = alpha_label;

                size++;
            }
        }

        delete [] variables;
    }

    delete e;
}


/**********************************************************************************************/
/* Performs alpha-expansion for  regular grid graph for case when energy terms are NOT        */
/* specified by a function */
void Expansion::set_up_expansion_energy_G_ARRAY(int size, LabelType alpha_label,Energy *e,
                                                     Energy::Var *variables )
{
    int i,nPix,pix,x,y;


    for ( i = size - 1; i >= 0; i-- )
    {
        pix = m_lookupPixVar[i];
        y = pix/m_width;
        x = pix - y*m_width;

        m_lookupPixVar[pix] = i;

        if ( x < m_width - 1 )
        {
            nPix = pix + 1;

            if ( m_labeling[nPix] != alpha_label )
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label),
                              m_smoothcost(alpha_label,m_labeling[nPix]),
                              m_smoothcost(m_labeling[pix],alpha_label),
                              m_smoothcost(m_labeling[pix],m_labeling[nPix]));
            else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
                                 m_smoothcost(m_labeling[pix],alpha_label));
        }   

        if ( y < m_height - 1 )
        {
            nPix = pix + m_width;

            if ( m_labeling[nPix] != alpha_label )
                e ->add_term2(variables[i],variables[m_lookupPixVar[nPix]],
                              m_smoothcost(alpha_label,alpha_label) ,
                              m_smoothcost(alpha_label,m_labeling[nPix]),
                              m_smoothcost(m_labeling[pix],alpha_label) ,
                              m_smoothcost(m_labeling[pix],m_labeling[nPix]) );
            else   e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
                                 m_smoothcost(m_labeling[pix],alpha_label));
        }   
        if ( x > 0 )
        {
            nPix = pix - 1;
    
            if ( m_labeling[nPix] == alpha_label )
               e ->add_term1(variables[i],m_smoothcost(alpha_label,m_labeling[nPix]),
                                 m_smoothcost(m_labeling[pix],alpha_label) );
        }   

        if ( y > 0 )
        {
            nPix = pix - m_width;
    
            if ( m_labeling[nPix] == alpha_label )
               e ->add_term1(variables[i],m_smoothcost(alpha_label,alpha_label),
                                 m_smoothcost(m_labeling[pix],alpha_label));
        }   
            
    }
    
}

/**************************************************************************************/

GraphCuts::EnergyType Expansion::oneExpansionIteration()
{
    int next;
   
    terminateOnError( m_dataType == NONE,"You have to set up the data cost before running optimization");
    terminateOnError( m_smoothType == NONE,"You have to set up the smoothness cost before running optimization");

    if (m_random_label_order) scramble_label_table();
    

    for (next = 0;  next < m_nLabels;  next++ )
        perform_alpha_expansion(m_labelTable[next]);
    
    
    return(dataEnergy()+smoothnessEnergy());
}


/**************************************************************************************/

void Expansion::optimizeAlg(int nIterations)
{
    expansion(nIterations);
}
