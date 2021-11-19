import os, sys
import argparse
import json
from Bio import SeqIO
import math
import csv

def log(num):
    if num == 0:
        return -math.inf
    else:
        return math.log(num)

def init(config):

    emit_p = {}
    start_p = {}
    trans_p = {}
    states = ['inter', 'start', 'mid', 'stop']

    start_p = {'inter':1, 'start':0, 'mid':0, 'stop':0}

    # set emission probabilities
    # ex: emit_p[state][observation]

    for state in states:
        emit_p[state] = config[state]

    # print(emit_p)

    # transition probabilites
    # ex: trans_p[from state][to state]

    inter_len = config['avg_len']['intergenic']
    gen_len = config['avg_len']['genic']/3

    # inter_len = 2
    # gen_len = 2

    trans_p['inter'] = {'inter': (inter_len-1)/inter_len, 'start':1/inter_len, 'mid':0, 'stop':0}
    trans_p['start'] = {'inter':0, 'start':0, 'mid':1, 'stop':0}
    trans_p['mid'] = {'inter':0, 'start':0, 'mid':(gen_len-1)/gen_len, 'stop':1/gen_len}
    trans_p['stop'] = {'inter':1, 'start':0, 'mid':0, 'stop':0}

    return states, start_p, emit_p, trans_p

def viterbi(contig_name, contig_seq, states, start_p, emit_p, trans_p):
    # store probabilities and pointers for each step of viterbi
    V = [{}]
    
    obs = str(contig_seq)
    # print(obs)

    obs = "GGATGAAATAA"


    # initiallize first observation at time 0
    for state in states:
        if state == 'inter':
            V[0][state] = {"prob": log(start_p[state]) + log(emit_p[state][obs[0]]), "prev": None, "loc": 0, "obs":obs[0], "is_last":False}
        else:
            V[0][state] = {"prob": -math.inf, "prev":None, "loc":2, "obs":obs[0:3], "is_last":False}

    # run viterbi for the rest of the sequence
    for t in range(1, len(obs)):
        V.append({})

        for new_state in states:
            
            max_state = 'inter'
            prev_prob = V[t-1][max_state]["prob"]

            if new_state == "inter":
                new_loc = max_loc = V[t-1][max_state]["loc"] + 1

                # if you are at the end of the sequence
                if new_loc+1 > len(obs):
                    new_loc = max_obs = None
                    new_prob = max_prob = -math.inf
                    # print(prev_state, new_state, new_loc, new_obs, new_prob)
                else:
                    new_obs = max_obs = obs[max_loc]
                    max_prob = prev_prob + log(trans_p[max_state][new_state]) + log(emit_p[new_state][new_obs])
                    # print(prev_state, new_state, new_loc, new_obs, new_prob)
            # for non-inter states
            else: 
                new_loc = max_loc = V[t-1][max_state]["loc"] + 3
                # reach end of seq
                if new_loc+1 > len(obs):
                    new_obs = None
                    new_prob = -math.inf
                    #print(prev_state, new_state, new_loc, new_obs, new_prob)
                else:
                    new_obs = max_obs = obs[max_loc-2: max_loc+1]
                    max_prob = prev_prob + log(trans_p[max_state][new_state]) + log(emit_p[new_state][max_obs])
                    #print(prev_state, new_state, new_loc, new_obs, new_prob)
            #print(max_loc, max_obs, max_prob, max_state)

            # for prev_states other than inter
            for prev_state in states[1:]:
                
                prev_prob = V[t-1][prev_state]["prob"]
      
                # get new obs for each possible new state
                if new_state == "inter":
                    new_loc = V[t-1][prev_state]["loc"] + 1
                    # if out of bounds
                    if new_loc+1 > len(obs):
                        new_obs = None
                        new_prob = -math.inf
                        #print(prev_state, new_state, new_loc, new_obs, new_prob)
                        
                    else:
                        new_obs = obs[new_loc]
                        new_prob = prev_prob + log(trans_p[prev_state][new_state]) + log(emit_p[new_state][new_obs])
                        #print(prev_state, new_state, new_loc, new_obs, new_prob)
                else:
                    new_loc = V[t-1][prev_state]["loc"] + 3
                    # reach end of seq
                    if new_loc+1 > len(obs):
                        new_obs = None
                        new_prob = -math.inf
                        #print(prev_state, new_state, new_loc, new_obs, new_prob)
                    else:
                        new_obs = obs[new_loc-2: new_loc+1]
                        new_prob = prev_prob + log(trans_p[prev_state][new_state]) + log(emit_p[new_state][new_obs])
                        #print(prev_state, new_state, new_loc, new_obs, new_prob)

                # get maximum probability from prev timestep
                # for each prev_state
                if new_prob > max_prob:
                    max_loc = new_loc
                    max_obs = new_obs
                    max_state = prev_state
                    max_prob = new_prob
            

            # if this obs is the last one of the sequence
            if max_loc == len(obs)-1:
                V[t][new_state] = {"prob": max_prob, "prev": max_state, "loc": max_loc, "obs":max_obs, "is_last":True}
            else:
                V[t][new_state] = {"prob": max_prob, "prev": max_state, "loc": max_loc, "obs":max_obs, "is_last":False}
            # print(max_prob, max_loc, max_obs)

        # stop loop if all observations are null (all states are past end of sequence)
        all_obs = []
        for state, info in V[t].items():
            all_obs.append(info['obs'])
        if all_obs.count(None) == 4:
            V.pop(t)
            break

    return(contig_name, V)

def traceback(V):
    # traceback
    # find all probabilities for is_last and get max
    # winner = {}
    path = []
    max_prob = -math.inf
    max_state = None
    max_i = None

    for i in range(0, len(V)):
        for state, info in V[i].items():
            if info['is_last']:
                print('found end', info['loc'], info['prob'])
                if info['prob'] > max_prob:
                    max_prob = info["prob"]
                    max_state = state
                    max_i = i
    # print(max_state, max_prob, max_i)

    # path.append(max_state)

    prev_state = max_state
    # traceback from max
    while prev_state != None:
        path.insert(0,[prev_state, V[max_i][prev_state]['loc']])
        prev_state = V[max_i][prev_state]['prev']
        max_i -= 1

    for i in path:
        if i[0] == 'start':
            print(i)
    # if 'start' in path:
    #     indices = [i for i, x in enumerate(path) if x == "start"]
    #     print("there is a start", path.count('start'), indices)
    # print('inters:', path.count('inter'))
    return(path)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', help='add the fasta file for prediction')
    parser.add_argument('-c','--config', help='add the config file')
    parser.add_argument('-o', '--output', help='add output file')
    args = parser.parse_args()
    fasta_file = args.fasta
    config_file = args.config
    output_file = args.output

    fasta_dict = SeqIO.to_dict(SeqIO.parse(open(fasta_file), 'fasta'))


    config= {}
    with open(config_file, 'r') as f:
        config = json.load(f)

    # print(config.keys())

    # initialize states and start probabilities
    states, start_p, emit_p, trans_p = init(config)

    # print(trans_p)
    # print(json.dumps(trans_p, indent = 4))

   

    # for contig in fasta_dict:
    #     contig_name = fasta_dict[contig].id
    #     contig_seq = fasta_dict[contig].seq

    contig_name = fasta_dict['DN38.contig00002'].id
    contig_seq = fasta_dict['DN38.contig00002'].seq

    contig_name, V = viterbi(contig_name, contig_seq, states, start_p, emit_p, trans_p)
            
        # print(json.dumps(V, indent=4))
    # for line in V:
    #     print(json.dumps(line, indent=4))

    path = traceback(V)
    # print(path)

    # generate gff3 file
    gff3 = []
    for i in path:
        if i[0] == 'start':
            start = i[1] - 1
        if i[0] == 'stop':
            stop = i[1] + 1
            gff3.append([contig_name, 'ena', start, stop, '.', '+', '0', '.'])

    # print(gff3)

    with open(output_file, 'w', newline='') as file:
        csvwriter = csv.writer(file, delimiter='\t')
        for line in gff3:
            csvwriter.writerow(line)




if __name__ == "__main__":
    main()