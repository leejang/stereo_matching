/*
 * File: LoopyBP.cpp
 *
 * motivated from the original MRF 2.2 Library
 * hacked by Jangwon Lee, 11/2016
 * Email: leejang@indiana.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <limits.h>

#include "LoopyBP.h"

using namespace std;

#define get_data_cost(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define get_smoothness_cost(l1,l2)  m_V[(l1)*m_nLabels+(l2)]
#define TRUNCATE_MIN(a,b) { if ((a) > (b)) (a) = (b); }


// inline functions from MRF 2.2
// to copy and add values the array of cost values
inline void CopyVector(MRF::CostVal* to, MRF::CostVal* from, int K)
{
    MRF::CostVal* to_finish = to + K;
    do
	{
	    *to ++ = *from ++;
	} while (to < to_finish);
}

inline void AddVector(MRF::CostVal* to, MRF::CostVal* from, int K)
{
    MRF::CostVal* to_finish = to + K;
    do
	{
	    *to ++ += *from ++;
	} while (to < to_finish);
}
inline MRF::CostVal min(MRF::CostVal a, MRF::CostVal b) {
    return a < b ? a : b;
}
inline MRF::CostVal SubtractMin(MRF::CostVal *D, int K)
{
    int k;
    MRF::CostVal delta;
	
    delta = D[0];
    for (k=1; k<K; k++) TRUNCATE_MIN(delta, D[k]);
    for (k=0; k<K; k++) D[k] -= delta;

    return delta;
}

LoopyBP::LoopyBP(int width, int height, int nLabels,EnergyFunction *eng):MRF(width,height,nLabels,eng)
{
    m_D = NULL;
    m_V = NULL;

    m_messages = NULL;
    m_messageArraySizeInBytes = 0;

    m_answer = new Label[m_nPixels];
}

LoopyBP::~LoopyBP()
{
    delete [] m_answer;
    delete [] m_V;

    if (m_messages)
        delete [] m_messages;
}

void LoopyBP::clearAnswer()
{
    // initialize all labels
    memset(m_answer, 0, m_nPixels*sizeof(Label));

    // initialize all messages
    if (m_messages) {
         memset(m_messages, 0, m_messageArraySizeInBytes);
    }
}

MRF::EnergyVal LoopyBP::dataEnergy()
{
    EnergyVal eng = (EnergyVal)0;

    for (int i = 0; i < m_nPixels; i++) {
        eng = eng + get_data_cost(i, m_answer[i]);
    }
    return (eng);
}

MRF::EnergyVal LoopyBP::smoothnessEnergy()
{
    EnergyVal eng = (EnergyVal)0;
    // pixel number
    unsigned int pix = 0;

    if (m_smoothType != FUNCTION ) {

        // x - direction
        for (int y = 0; y < m_height; y++) {
            for (int x = 1; x < m_width; x++ ) {
                pix = x + y * m_width;
                eng = eng + get_smoothness_cost(m_answer[pix], m_answer[pix-1]);
            }
        }
        // y - direction
        for (int y = 1; y < m_height; y++) {
            for (int x = 0; x < m_width; x++ ) {
                pix = x + y * m_width;
                eng = eng + get_smoothness_cost(m_answer[pix], m_answer[pix - m_width]);
            }
        }
    } else {
        // function type
    }

    return (eng);
}

void LoopyBP::setData(CostVal* data)
{
    m_D = data;
}

void LoopyBP::setSmoothness(SmoothCostGeneralFn cost)
{
}

// Truncated L1
void LoopyBP::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
{
    assert(lambda > 0);

    int ki, kj;
    CostVal cost;

    m_V = new CostVal[m_nLabels*m_nLabels];

    for (ki = 0; ki < m_nLabels; ki++) {
        for (kj = ki; kj < m_nLabels; kj++) {
                // L1
		cost = (CostVal) (kj - ki);
		if (cost > smoothMax)
                    cost = smoothMax;

		m_V[ki*m_nLabels + kj] = m_V[kj*m_nLabels + ki] = cost*lambda;
	    }
    }

    m_smoothMax = smoothMax;
    m_lambda = lambda;

}

#if 0
// Pott Model
void LoopyBP::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
{
    assert(lambda > 0);

    int ki, kj;
    CostVal cost;

    m_V = new CostVal[m_nLabels*m_nLabels];

    for (ki = 0; ki < m_nLabels; ki++) {
        for (kj = ki; kj < m_nLabels; kj++) {

                if (kj == ki)
		  cost = (CostVal)0;
                else
		  cost = (CostVal)10;

		m_V[ki*m_nLabels + kj] = m_V[kj*m_nLabels + ki] = cost;
	    }
    }

    m_smoothMax = smoothMax;
    m_lambda = lambda;

}
#endif

void LoopyBP::initializeAlg()
{
    int messageNum = 4 * m_nPixels * m_nLabels;
    m_messageArraySizeInBytes = messageNum*sizeof(CostVal);

    m_messages = new CostVal[messageNum];
    memset(m_messages, 0, messageNum*sizeof(CostVal));
}

//first implementation
void LoopyBP::sendMsg(int x, int y, MRF::CostVal* data_i, MRF::CostVal* msg_i, MRF::CostVal *msgs, DIRECTION direction)
{
    CostVal new_msg[m_nLabels];

    CostVal min_cost;

    for (int i = 0; i < m_nLabels; i++) {
        // maximum value of unsigned int
        min_cost = INT_MAX;

        for (int j = 0; j < m_nLabels; j++) {

            CopyVector(msg_i, data_i + i + j, 1);
            //cout << "1 " << endl;
            // Exclude the incoming message direction
            if ((direction != RIGHT) && (x > 0)) {
                AddVector(msg_i, msgs - 2*m_nLabels + i + j, 1); 
            }
            if ((direction != LEFT) && (x < m_width - 1)) {
                AddVector(msg_i, msgs + i + j, 1);
            }
            if ((direction != UP) && (y > 0)) {
                AddVector(msg_i, msgs - (2*m_width - 1)*m_nLabels + i + j, 1); 
            }
            if ((direction != DOWN) && (y < m_height - 1)) {
                AddVector(msg_i, msgs + m_nLabels+ i + j, 1);
            }

            if (msg_i[0] < 0)
                msg_i[0] = min_cost;

            if (min_cost > msg_i[0])
                min_cost = msg_i[0];
        }
        new_msg[i] = min_cost;
    }

    // update messages
    for (int i = 0; i < m_nLabels; i++) {
        //cout << "i: " << i << ",  new_msg[i]: " << new_msg[i] << endl;
        switch(direction) {
            case RIGHT:
              if (x > 0) msgs[i - 2*m_nLabels] = new_msg[i];
            break;
            case LEFT:
              if (x < m_width - 1) msgs[i] = new_msg[i];
            break;
            case UP:
              if (y > 0) msgs[i - (2*m_width - 1)*m_nLabels] = new_msg[i];
            break;
            case DOWN:
              if (y < m_height - 1) msgs[i + m_nLabels] = new_msg[i]; 
            break;
        }
    }

}

// modified version to speed up
void LoopyBP::sendMsg(MRF::CostVal* msg_i, MRF::CostVal *msgs, DIRECTION direction)
{
    switch(direction) {
        case RIGHT:
            AddVector(msg_i, msgs - 2*m_nLabels, m_nLabels);
        break;
        case LEFT:
            AddVector(msg_i, msgs, m_nLabels);
        break;
        case UP:
            AddVector(msg_i, msgs - (2*m_width - 1)*m_nLabels, m_nLabels);
        break;
        case DOWN:
            AddVector(msg_i, msgs + m_nLabels, m_nLabels);
        break;
    }

}

// update Messages
void LoopyBP::updateMsgs(MRF::CostVal* msgs, MRF::CostVal* msg_i)
{
    CostVal min_cost;

    int idx = 0;
    min_cost = msgs[0] = msg_i[0] - msgs[0];
    for (idx = 1; idx < m_nLabels; idx++) {
        msgs[idx] = msg_i[idx] - msgs[idx];
        TRUNCATE_MIN(min_cost, msgs[idx]);
        TRUNCATE_MIN(msgs[idx], msgs[idx - 1] + m_lambda);
    }

    msgs[--idx] -= min_cost;
    TRUNCATE_MIN(msgs[idx], m_lambda*m_smoothMax);

    for (idx--; idx >= 0; idx--) {
        msgs[idx] -= min_cost;
        TRUNCATE_MIN(msgs[idx], msgs[idx + 1] + m_lambda);
        TRUNCATE_MIN(msgs[idx], m_lambda*m_smoothMax);
    }
    //cout << "update Msgs" << endl;
}

void LoopyBP::optimizeAlg(int nIterations)
{
    CostVal msg_i[m_nLabels];

    // pointers indicate data and messages
    // data: data_cost + smoothness_cost
    CostVal* data_pt;
    CostVal* msgs_pt;

    int idx = 0;

    while (nIterations > 0) {

        idx = 0;
        data_pt = m_D;
        msgs_pt = m_messages;

        for (int y = 0; y < m_height; y++) {
            for (int x = 0; x < m_width; x++) {
                CopyVector(msg_i, data_pt, m_nLabels);
                if (x > 0) sendMsg(msg_i, msgs_pt, RIGHT);
                if (x < m_width - 1) sendMsg(msg_i, msgs_pt, LEFT);
                if (y > 0) sendMsg(msg_i, msgs_pt, UP);
                if (y < m_height - 1) sendMsg(msg_i, msgs_pt, DOWN);

                // update messages
                if (x < m_width - 1) {
                    updateMsgs(msgs_pt, msg_i);
                }
                if (y < m_height - 1) {
                    updateMsgs(msgs_pt + m_nLabels, msg_i);
                }

                data_pt += m_nLabels;
                msgs_pt += 2*m_nLabels;
                idx++;
            }
        }
 
        idx --;
        data_pt -= m_nLabels;
        msgs_pt -= 2*m_nLabels;
    
        for (int y = m_height - 1; y >= 0; y--) {
            for (int x = m_width - 1; x>= 0; x--) {

                CopyVector(msg_i, data_pt, m_nLabels);
                if (x > 0) sendMsg(msg_i, msgs_pt, RIGHT);
                if (x < m_width - 1) sendMsg(msg_i, msgs_pt, LEFT);
                if (y > 0) sendMsg(msg_i, msgs_pt, UP);
                if (y < m_height - 1) sendMsg(msg_i, msgs_pt, DOWN);

                SubtractMin(msg_i, m_nLabels);

                // update messages
                if (x > 0) {
                    updateMsgs(msgs_pt - 2*m_nLabels, msg_i);
                }
                if (y > 0) {
                    updateMsgs(msgs_pt - (2*m_width - 1)*m_nLabels, msg_i);
                }

                data_pt -= m_nLabels;
                msgs_pt -= 2*m_nLabels;
                idx--;
            }
        }
        nIterations--;
    } // main while loop

   // get answers
   getAnswers();
}

void LoopyBP::getAnswers()
{
    CostVal* data_i;
    data_i = new CostVal[m_nLabels];

    // pointers indicate data and messages
    CostVal* data_pt;
    CostVal* message_pt;

    // data and messages
    data_pt = m_D;
    message_pt = m_messages;

    // index
    unsigned int idx = 0;
    // label index
    unsigned int kj = 0;
    
    CostVal min_cost;

    for (int y = 0; y < m_height; y++) {	    
        for (int x = 0; x < m_width; x++) {
            // copy the data_i to handle
            CopyVector(data_i, data_pt, m_nLabels);
            // Right
            if (x > 0) {
                kj = m_answer[idx - 1];
                for (int label = 0; label < m_nLabels; label++) {
                    data_i[label] += m_V[kj * m_nLabels + label]; 
                }
            }
            // Up
            if (y > 0) {
                kj = m_answer[idx - m_width];
                for (int label = 0; label < m_nLabels; label++) {
                    data_i[label] += m_V[kj * m_nLabels + label]; 
                }
            }
            // Left
            if (x < m_width - 1) {
                AddVector(data_i, message_pt, m_nLabels);
            }
            // Down
            if (y < m_height - 1) {
                AddVector(data_i, message_pt + m_nLabels, m_nLabels);
            }

            min_cost = data_i[0];
            m_answer[idx] = 0;

            for (int label = 1; label < m_nLabels; label++) {
                // update when we have min cost
                if (min_cost > data_i[label]) {
                    min_cost = data_i[label];
                    // assign new label
                    m_answer[idx] = label;
                }
            }

            data_pt += m_nLabels;
            message_pt += 2*m_nLabels;
            idx++;
	}
    }

   delete [] data_i;
}
