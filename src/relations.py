# This is the file containing the functions handling the relations list

from utils import invmod
from polynomial_functions import *

def find_cycle_fp(graph, init):
    if init[0] == 1: return DFS_fp(graph, init)

    path_p1_to_1 = DFS_fp(graph, [1, init[0]])

    if path_p1_to_1 is not None:
        path_p2_to_1 = DFS_fp(graph, [1, init[1]])

        if path_p2_to_1 is not None:
            index1, index2 = len(path_p1_to_1)-1, len(path_p2_to_1)-1
            while path_p1_to_1[index1] == path_p2_to_1[index2]:
                index1 -= 1
                index2 -= 1
            return path_p2_to_1[:index2+2] + path_p1_to_1[:index1+1][::-1]

    return DFS_fp(graph, init)

# Finds cycles, aka relations, in the graph of relations with one or two large primes in the rational element
def DFS_fp(graph, init):
    path = [init[1]]
    stack = []
    neighbors = [init[1]]
    
    # Initialize the potential paths
    tmp = graph[init[1]]

    if init[0] in tmp:
        return path + [init[0]]
    
    neighbors += tmp
    
    stack.append(neighbors)
    
    # While there is a potential path that has not been explored
    while len(stack):
        while len(stack) and len(stack[-1]) == 1:
            path.pop()
            stack.pop()
        if not len(stack): break

        next_node = stack[-1][-1]
        stack[-1].pop()
        path.append(next_node)
        neighbors = [next_node]

        if init[0] in graph[next_node]:
            return path + [init[0]]
        
        parent = stack[-1][0]

        for i in range(len(graph[next_node])):
            if graph[next_node][i] != parent: neighbors.append(graph[next_node][i])
            else:
                neighbors += graph[next_node][i+1:]
                break

        stack.append(neighbors)
        
# Given a cycle in the fp graph, constructs the full relation
# The graph is always acyclic: when we find a new cycle, the partial relation we just found is not added
# Thus, one new partial relation can only create one cycle at most
def combine_fp(pair, pairs_used, fp, path, divide_leading, g):
    tmp = [i for i in pair]
    for i in range(len(path)-1):
        firstp, secondp = min(path[i],path[i+1]), max(path[i],path[i+1])

        combine = fp[firstp][secondp]
        tmp2 = [div_poly(poly_prod(tmp[0], combine[0]), g),
                tmp[1]*combine[1],
                [[1, 1]],
                tmp[3]*combine[3],
                tmp[4] + [path[i]],
                poly_prod(tmp[5],combine[5]),
                tmp[6]+1]
        for p in range(len(divide_leading)):
            tmp2.append(tmp[7+p] or combine[7+p])

        tmp = [i for i in tmp2]
    pairs_used.append(tmp)
    
# Master function, keeps the partial_relations and possible_smooth lists sorted, find the matching large primes, create the full relations from large primes
def handle_large_fp(pair, large_primes, pairs_used, fp, graph_fp, size_fp, parent_fp, g, divide_leading, cycle_len, full_found, partial_found_fp):
    pair[4] = [large_primes[3]]

    if large_primes[3] == large_primes[4]:
        pairs_used.append(pair)
        full_found += 1

    elif not bool(fp):
        pair[4].append(large_primes[4])

        small_p, big_p = large_primes[3], large_primes[4]

        fp[small_p] = {}
        fp[small_p][big_p] = pair

        graph_fp[small_p] = [big_p]
        graph_fp[big_p] = [small_p]

        parent_fp[small_p] = small_p
        parent_fp[big_p] = small_p

        size_fp += 1

    else:
        small_p, big_p = large_primes[3], large_primes[4]

        flag_small_prime = small_p in graph_fp
        
        # If the smallest prime has not already been seen
        if not flag_small_prime:
            graph_fp[small_p] = [big_p]

            # Find where in partial relations it is needed to insert the new partial relation
            pair[4].append(big_p)

            fp[small_p] = {}
            fp[small_p][big_p] = pair

            size_fp += 1
            
            # Find if largest prime has already been seen
            flag_big_prime = big_p in graph_fp

            if flag_big_prime: # if seen add the smallest prime to the list of links of the largest prime
                graph_fp[big_p].append(small_p)

                parent_fp[small_p] = parent_fp[big_p]

            else: # If it has not been seen, append it
                graph_fp[big_p] = [small_p]

                parent_fp[small_p] = small_p
                parent_fp[big_p] = small_p
            
        # If the smallest prime has already been seen
        else:
            flag_big_prime = big_p in graph_fp
            
            # If big prime has not been seen
            if not flag_big_prime:
                graph_fp[small_p].append(big_p) # Append the largest prime to the list of links of the smallest prime
                graph_fp[big_p] = [small_p] # Create the large prime in the graph

                parent_fp[big_p] = parent_fp[small_p]

                pair[4].append(big_p)
                # Find where to insert the partial relation
                if small_p not in fp:
                    fp[small_p] = {}

                fp[small_p][big_p] = pair

                size_fp += 1

            # If the largest prime has been seen, ie if both primes have already been seen
            else:

                parent_small_p = parent_fp[small_p]
                while parent_fp[parent_small_p] != parent_small_p:
                    parent_fp[parent_small_p], parent_small_p = parent_fp[parent_fp[parent_small_p]], parent_fp[parent_small_p]

                parent_big_p = parent_fp[big_p]
                while parent_fp[parent_big_p] != parent_big_p:
                    parent_fp[parent_big_p], parent_big_p = parent_fp[parent_fp[parent_big_p]], parent_fp[parent_big_p]

                if parent_small_p != parent_big_p:
                    graph_fp[small_p].append(big_p)
                    graph_fp[big_p].append(small_p)

                    parent_fp[parent_big_p] = parent_small_p

                    pair[4].append(big_p)
                    # Find where to insert the partial relation
                    if small_p not in fp:
                        fp[small_p] = {}
                        
                    fp[small_p][big_p] = pair

                    size_fp += 1
                else:
                    path_cycle = find_cycle_fp(graph_fp, [small_p, big_p])
                    if len(path_cycle) < 11: cycle_len[len(path_cycle)-2] += 1
                    else: cycle_len[-1] += 1
                    combine_fp(pair, pairs_used, fp, path_cycle, divide_leading, g)
                    partial_found_fp += 1

    return pairs_used, fp, graph_fp, size_fp, parent_fp, cycle_len, full_found, partial_found_fp

def find_cycle_pf(graph, init):
    if init[0] == (1, 1): return DFS_pf(graph, init)

    path_p1_to_1 = DFS_pf(graph, [(1, 1), init[0]])

    if path_p1_to_1 is not None:
        path_p2_to_1 = DFS_pf(graph, [(1, 1), init[1]])

        if path_p2_to_1 is not None:
            index1, index2 = len(path_p1_to_1)-1, len(path_p2_to_1)-1
            while path_p1_to_1[index1] == path_p2_to_1[index2]:
                index1 -= 1
                index2 -= 1
            return path_p2_to_1[:index2+2] + path_p1_to_1[:index1+1][::-1]

    return DFS_pf(graph, init)


# Finds cycles, aka relations, in the graph of relations with one or two large primes if the norm of the algebraic element
def DFS_pf(graph, init):
    path = [init[1]]
    stack = []
    neighbors = [init[1]]
    
    tmp = graph[init[1]]

    if init[0][0] in tmp and init[0] in tmp[init[0][0]]:
        return path + [init[0]]

    # Initialize the potential paths
    for key in tmp.keys():
        neighbors += tmp[key]
    
    stack.append(neighbors)
    
    # While there is a potential path that has not been explored
    while len(stack):
        while len(stack) and len(stack[-1]) == 1:
            path.pop()
            stack.pop()
        if not len(stack): break
        next_node = stack[-1][-1]
        stack[-1].pop()
        path.append(next_node)
        neighbors = [next_node]

        tmp = graph[next_node]

        if init[0][0] in tmp and init[0] in tmp[init[0][0]]:
            return path + [init[0]]
        
        parent = stack[-1][0]

        for key in tmp.keys():
            if key == parent[0]:
                for i in range(len(tmp[key])):
                    if tmp[key][i] != parent: neighbors.append(tmp[key][i])
                    else:
                        neighbors += tmp[key][i+1:]
                        break

            else:
                neighbors += tmp[key]

        stack.append(neighbors)


# Given a cycle in the pf graph, constructs the full relation
# The graph is always acyclic: when we find a new cycle, the partial relation we just found is not added
# Thus, one new partial relation can only create one cycle at most
def combine_pf(pair, pairs_used, pf, path, divide_leading, g):
    tmp = [i for i in pair]
    for i in range(len(path)-1):
        if path[i][0] <= path[i+1][0]: firstp, secondp = path[i], path[i+1]
        else: firstp, secondp = path[i+1], path[i]

        combine = pf[firstp][secondp]
        tmp2 = [div_poly(poly_prod(tmp[0], combine[0]), g),
                tmp[1]*combine[1],
                tmp[2] + [path[i]],
                tmp[3]*combine[3],
                [1],
                poly_prod(tmp[5], combine[5]),
                tmp[6]+1]
        for p in range(len(divide_leading)):
            tmp2.append(tmp[7+p] or combine[7+p])

        tmp = [i for i in tmp2]

    pairs_used.append(tmp)

def handle_large_pf(pair, large_primes, pairs_used, pf, graph_pf, size_pf, parent_pf, g, divide_leading, cycle_len, full_found, partial_found_pf):
    if large_primes[1] == 1:
        tmp = (1, 1) # [smallest_prime, r_smallest]
        pair[2] = [tmp]

    else:
        tmp = (large_primes[1], pair[5][1]*invmod(-pair[5][0], large_primes[1])%large_primes[1]) # [smallest_prime, r_smallest]
        pair[2] = [tmp]

    tmp2 = (large_primes[2], pair[5][1]*invmod(-pair[5][0], large_primes[2])%large_primes[2]) # [largest_prime, r_largest]

    if large_primes[1] == large_primes[2]:
        pairs_used.append(pair)
        full_found += 1

    elif not len(pf):
        pair[2].append(tmp2)
        pf[tmp] = {}
        pf[tmp][tmp2] = pair

        graph_pf[tmp] = {}
        graph_pf[tmp][tmp2[0]] = [tmp2]
        graph_pf[tmp2] = {}
        graph_pf[tmp2][tmp[0]] = [tmp]

        parent_pf[tmp] = tmp
        parent_pf[tmp2] = tmp

    else:
        flag_small_prime = tmp in graph_pf
        
        # If [smallest_prime, r_smallest] has not already been seen
        if not flag_small_prime:
            graph_pf[tmp] = {}
            graph_pf[tmp][tmp2[0]] = [tmp2] # Create smallest prime in graph

            pair[2].append(tmp2)
            pf[tmp] = {}
            pf[tmp][tmp2] = pair
            
            # Find if largest prime has already been seen
            flag_big_prime = tmp2 in graph_pf

            if flag_big_prime:
                if tmp[0] in graph_pf[tmp2]: # If seen, add the smallest prime to the list of links of the largest prime
                    graph_pf[tmp2][tmp[0]].append(tmp)
                else:
                    graph_pf[tmp2][tmp[0]] = [tmp]

                parent_pf[tmp] = tmp2

            else: # If it has not been seen, create it in the graph
                graph_pf[tmp2] = {}
                graph_pf[tmp2][tmp[0]] = [tmp]

                parent_pf[tmp] = tmp
                parent_pf[tmp2] = tmp

        # If [smallest_prime, r_smallest] has already be seen
        # We look for [largest_prime, r_largest]
        else:
            flag_big_prime = tmp2 in graph_pf

            if not flag_big_prime: # [largest_prime, r_biggest] not in graph
                if tmp2[0] in graph_pf[tmp]:
                    graph_pf[tmp][tmp2[0]].append(tmp2) # Create the link from small to big prime
                else:
                    graph_pf[tmp][tmp2[0]] = [tmp2] # Create the link from small to big prime

                graph_pf[tmp2] = {}
                graph_pf[tmp2][tmp[0]] = [tmp] # Create the large prime in the graph

                parent_pf[tmp2] = tmp

                pair[2].append(tmp2)
                # Find where to insert the partial relation
                if tmp not in pf:
                    pf[tmp] = {}

                pf[tmp][tmp2] = pair

                size_pf += 1

            else:
                parent_tmp = parent_pf[tmp]
                while parent_pf[parent_tmp] != parent_tmp:
                    parent_pf[parent_tmp], parent_tmp = parent_pf[parent_pf[parent_tmp]], parent_pf[parent_tmp]

                parent_tmp2 = parent_pf[tmp2]
                while parent_pf[parent_tmp2] != parent_tmp2:
                    parent_pf[parent_tmp2], parent_tmp2 = parent_pf[parent_pf[parent_tmp2]], parent_pf[parent_tmp2]

                if parent_tmp != parent_tmp2:
                    if tmp2[0] in graph_pf[tmp]:
                        graph_pf[tmp][tmp2[0]].append(tmp2)
                    else:
                        graph_pf[tmp][tmp2[0]] = [tmp2]

                    if tmp[0] in graph_pf[tmp2]:
                        graph_pf[tmp2][tmp[0]].append(tmp)
                    else:
                        graph_pf[tmp2][tmp[0]] = [tmp]

                    parent_pf[parent_tmp2] = parent_tmp

                    pair[2].append(tmp2)

                    if tmp not in pf:
                        pf[tmp] = {}
                        
                    pf[tmp][tmp2] = pair

                    size_pf += 1

                else:
                    path_cycle = find_cycle_pf(graph_pf, [tmp, tmp2])
                    if len(path_cycle) < 11: cycle_len[len(path_cycle)-2] += 1
                    else: cycle_len[-1] += 1
                    combine_pf(pair, pairs_used, pf, path_cycle, divide_leading, g)
                    partial_found_pf += 1

    return pairs_used, pf, graph_pf, size_pf, parent_pf, cycle_len, full_found, partial_found_pf