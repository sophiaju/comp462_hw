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
    V = []
    V.append({})
    
    obs = str(contig_seq)
    # print(obs)

    # obs = "GGATGATGGGGTAAC"


    # initiallize first observation at time 0
    for state in states:
        if state == 'inter':
            V[0][state] = {"prob": log(start_p[state]) + log(emit_p[state][obs[0]]), "prev": None, "prev_t": None, "loc": 0, "obs":obs[0]}
        else:
            V[0][state] = {"prob": -math.inf, "prev":None, "prev_t": None, "loc":2, "obs":obs[0:3]}

    # run viterbi for the rest of the sequence
    for t in range(1, len(obs)):
        V.append({})
        # if t < 5:
        #     for line in V:
        #         print(json.dumps(line, indent=4))

        for new_state in states:
            # print(new_state)
            
            # check if non-inter state loc and intergenic loc are the same, if not then make copy (wait)
            if new_state == 'start':
                if math.isfinite(V[t-1]['start']['prob']) and (V[t-1]['start']['loc'] != V[t-1]['inter']['loc']):
                    # print(t, V[t-1]['start']['prob'], 'start')
                    V[t]['start'] = V[t-1]['start']
                    V[t][new_state]['note'] = 'copied last state, curr state: ' + str(t)
                    V[t]['mid'] = V[t-1]['mid']
                    V[t]['mid']['note'] = 'copied last state, curr state: ' + str(t)
                    V[t]['stop'] = V[t-1]['stop']
                    V[t]['stop']['note'] = 'copied last state, curr state: ' + str(t)
                    break
            if new_state == 'mid' and (math.isfinite(V[t-1]['mid']['prob'] or math.isfinite(V[t-1]['start']['prob']))):
                # print(t, V[t-1]['mid']['prob'], 'mid')
                if (V[t-1]['start']['loc'] != V[t-1]['inter']['loc']):
                    # print((V[t-1]['start']['loc'], V[t-1]['inter']['loc']), 'mid')
                    V[t][new_state] = V[t-1][new_state]
                    V[t]['mid']['note'] = 'copied last state, curr state: ' + str(t)
                    V[t]['start'] = V[t-1]['start']
                    V[t]['start']['note'] = 'copied last state, curr state: ' + str(t)
                    V[t]['stop'] = V[t-1]['stop']
                    V[t]['stop']['note'] = 'copied last state, curr state: ' + str(t)
                    break
            if new_state == 'stop' and math.isfinite(V[t-1]['stop']['prob']):
                # print(t, V[t-1]['mid']['prob'])
                if (V[t-1]['mid']['loc'] != V[t-1]['inter']['loc']):
                    V[t][new_state] = V[t-1][new_state]
                    V[t][new_state]['note'] = 'copied last state, curr state: ' + str(t)
                    V[t]['start'] = V[t-1]['start']
                    V[t]['start']['note'] = 'copied last state, curr state: ' + str(t)
                    V[t]['mid'] = V[t-1]['mid']
                    V[t]['mid']['note'] = 'copied last state, curr state: ' + str(t)
                    break
            

            # if they are the same, do trans prob comparison
            max_prob_t = V[t-1]['inter']['prob'] + log(trans_p['inter'][new_state])
            max_state = 'inter'

            for prev_state in states[1:]:

                if new_state == 'inter' and prev_state == 'stop' and math.isfinite(V[t-1]['stop']['prob']) and (V[t-1]['stop']['loc'] != V[t-1]['inter']['loc']):
                    # print('dont compare stop to inter', V[t-1]['stop']['loc'], V[t-1]['inter']['loc'])
                    new_prob_t = -math.inf
                   
                else:
                    new_prob_t = V[t-1][prev_state]['prob'] + log(trans_p[prev_state][new_state])

                    # if new_state == 'inter':
                        # print(new_prob_t, prev_state,new_state, "new inter prob")
                
                if new_prob_t > max_prob_t:
                    max_prob_t = new_prob_t
                    max_state = prev_state

            prev_loc = V[t-1][max_state]['loc']

            # get new observation
            if new_state == 'inter':
                new_loc = prev_loc + 1
                new_obs = obs[new_loc]
            else:
                new_loc = prev_loc + 3
                new_obs = obs[new_loc-2: new_loc+1]
            
            # if end of sequence and there is wrong observation format
            if new_obs not in emit_p[new_state]:
                total_prob = -math.inf
            else:
                # get max probability with emit
                total_prob = max_prob_t + log(emit_p[new_state][new_obs])

            # update new column
            V[t][new_state] = {"prob": total_prob, "prev": max_state, "prev_t":t-1, "loc": new_loc, "obs":new_obs, 'note':'did not copy last state, curr state: ' + str(t)}
            # V[t][new_state]['note'] = 'did not copy last state, curr state: ' + str(t)

        # print(t)
        # print(json.dumps(V[t], indent=4))
    return(contig_name, V)

def traceback(V):
    # traceback
    # find all probabilities for is_last and get max
    # winner = {}
    path = []
    max_prob = -math.inf
    max_state = None
    prev_t = None
    prev = None
    max_loc = None
    for state, info in V[-1].items():
        if info['prob'] > max_prob:
            max_prob = info['prob']
            max_state = state
            prev_t = info['prev_t']
            prev = info['prev']
            max_loc = info['loc']
            
    path.append([max_state, max_loc])

    # traceback from max
    while prev:
        # print(prev, prev_t)
        # print(V[prev_t][prev])
        path.insert(0,[prev, V[prev_t][prev]['loc']])
        new_prev = V[prev_t][prev]['prev']
        new_prev_t = V[prev_t][prev]['prev_t']
        prev = new_prev
        prev_t = new_prev_t

    # for i in path:
    #     if i[0] == 'start':
    #         print(i)
    # if 'start' in path:
    #     indices = [i for i, x in enumerate(path) if x == "start"]
    #     print("there is a start", path.count('start'), indices)
    # print('inters:', path.count('inter'))
    return(path)

def gen_gff3(contig_name, path):
    gff3 = []
    for i in path:
        if i[0] == 'start':
            start = i[1] - 1
        if i[0] == 'stop':
            stop = i[1] + 1
            gff3.append([contig_name, 'ena', 'CDS', start, stop, '.', '+', '0', '.'])

    return(gff3)


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

    all_gff3 = []

    for contig in fasta_dict:
        contig_name = fasta_dict[contig].id
        contig_seq = fasta_dict[contig].seq

    # contig_name = fasta_dict['DN38.contig00002'].id
    # contig_seq = fasta_dict['DN38.contig00002'].seq

        contig_name, V = viterbi(contig_name, contig_seq, states, start_p, emit_p, trans_p)
                
        # print(V[0])
        # for line in V:
        #     print(json.dumps(line, indent=4))

        path = traceback(V)
        # print(path)

        # generate gff3 file

        all_gff3.append(gen_gff3(contig_name, path))

    with open(output_file, 'w', newline='') as file:
        csvwriter = csv.writer(file, delimiter='\t')
        for contig in all_gff3:
            for line in contig:
                csvwriter.writerow(line)




if __name__ == "__main__":
    main()