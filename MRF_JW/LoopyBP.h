/*
 * File: LoopyBP.h
 *
 * motivated from the original MRF 2.2 Library
 * hacked by Jangwon Lee, 11/2016
 * Email: leejang@indiana.edu
 */

#ifndef __LoopyBP_H__
#define __LoopyBP_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "mrf.h"

class LoopyBP : public MRF {
 public:
    LoopyBP(int width, int height, int nLabels, EnergyFunction *eng);
    ~LoopyBP();
    void clearAnswer();
    // Energy Functions
    EnergyVal dataEnergy();
    EnergyVal smoothnessEnergy();

    void setLabel(int pixel,Label label){m_answer[pixel] = label;};

    // Below functions are NOT used in this class
    void setNeighbors(int /*pix1*/, int /*pix2*/, CostVal /*weight*/){printf("Not implemented"); exit(1);}
    void setParameters(int /*numParam*/, void * /*param*/){printf("No optional parameters to set"); exit(1);}

    Label getLabel(int pixel){return(m_answer[pixel]);};
    Label* getAnswerPtr(){return(m_answer);};

 protected:
    // for Data cost function
    void setData(DataCostFn /*dcost*/){printf("Not implemented"); exit(1);}
    void setData(CostVal* data);
    // for Smoothness function
    // to assign specific weights for each node
    void setCues(CostVal* /*hCue*/, CostVal* /*vCue*/){printf("Not implemented"); exit(1);}
    void setSmoothness(CostVal* /*V*/){printf("Not implemented"); exit(1);}
    void setSmoothness(SmoothCostGeneralFn cost);
    void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda);
    // to run algorithm
    enum DIRECTION {RIGHT, LEFT, UP, DOWN};

    void initializeAlg();
    void sendMsg(int x, int y, MRF::CostVal* data_i, MRF::CostVal* msg_i, MRF::CostVal* msgs, DIRECTION direction);
    void sendMsg(MRF::CostVal* msg_i, MRF::CostVal* msgs, DIRECTION direction);
    void updateMsgs(MRF::CostVal* msgs, MRF::CostVal* msg_i);

    void optimizeAlg(int nIterations);
    // to get final answers
    void getAnswers();

 private:
    // stored data cost values
    CostVal *m_D;
    // stored smoothness cost values
    CostVal *m_V;

    // variables for smoothness
    CostVal m_smoothMax;
    CostVal m_lambda;

    // for answerr ans messages
    Label *m_answer;
    CostVal *m_messages;
    // size of messages
    int	m_messageArraySizeInBytes;
};

#endif /*  __LoopyBP_H__ */
