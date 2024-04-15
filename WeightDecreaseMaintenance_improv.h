#pragma once

#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

void WeightDecreaseMaintenance_improv_step1(int v1, int v2, weightTYPE w_new, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR, std::vector<affected_label>* CL,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) 
{

	for (int sl = 0; sl < 2; sl++) {
		if (sl == 1) {
			swap(v1, v2);
		}
		//第一次遍历L(a) 第二次遍历L(b)
		for (auto it : (*L)[v1]) {
			//即r(v)>=r(b) 从小到大排序
			if (it.vertex <= v2) {
				results_dynamic.emplace_back(pool_dynamic.enqueue([it, v2, L, PPR, w_new, CL] {
					
					//计算Query距离
					auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, it.vertex, v2); // query_result is {distance, common hub}
					if (query_result.first > it.distance + w_new) {
						mtx_595_1.lock();
						CL->push_back(affected_label{ v2 , it.vertex, it.distance + w_new });
						mtx_595_1.unlock();
					}
					else {
						//L(b)[v]
						auto search_result = search_sorted_two_hop_label((*L)[v2], it.vertex);
						if (search_result > it.distance + w_new && search_result != MAX_VALUE) {
							mtx_595_1.lock();
							CL->push_back(affected_label{ v2, it.vertex, it.distance + w_new });
							mtx_595_1.unlock();
						}
						if (query_result.second != it.vertex) {
							mtx_5952[v2].lock();
							PPR_insert(*PPR, v2, query_result.second, it.vertex);
							mtx_5952[v2].unlock();
						}
						if (query_result.second != v2) {
							mtx_5952[it.vertex].lock();
							PPR_insert(*PPR, it.vertex, query_result.second, v2);
							mtx_5952[it.vertex].unlock();
						}
					}

					return 1; }));
			}
		}
	}

	for (auto&& result : results_dynamic) {
		result.get();
	}
	std::vector<std::future<int>>().swap(results_dynamic);
}




void DIFFUSE(graph_v_of_v_idealID& instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR, std::vector<affected_label>& CL,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic)
{
    for (const auto& al : CL) 
	{
        results_dynamic.emplace_back(pool_dynamic.enqueue([L, PPR, al, &instance_graph] 
		{

            std::vector<weightTYPE> Dis(instance_graph.get_num_vertices(), -1);
            Dis[al.first] = al.dis;
            std::priority_queue<std::pair<weightTYPE, int>> Q;
            Q.emplace(al.dis, al.first);

            while (!Q.empty()) 
			{
                int x = Q.top().second;
                weightTYPE d_x = Q.top().first;
                Q.pop();

            
                auto& label_list = (*L)[x];
                auto it = std::lower_bound(label_list.begin(), label_list.end(), al.second, [](const two_hop_label_v1& l, int v) {
                    return l.vertex < v;
                });
                if (it != label_list.end() && it->vertex == al.second) {
                    it->distance = d_x;
                }
                else {
                    label_list.emplace(it, al.second, d_x);
                }

               
                for (int x_n : instance_graph.get_neighbors(x))
				{
                    if () 
					{
                        if (Dis[x_n] == -1) {
                            Dis[x_n] = instance_graph.graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L,x, x_n) + d_x;
                            Q.emplace(Dis[x_n], x_n);
                        }
                        else if (Dis[x_n] > instance_graph.graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L,x, x_n)(x, x_n) + d_x) {
                            Dis[x_n] = instance_graph.graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L,x, x_n)(x, x_n) + d_x;
                            Q.emplace(Dis[x_n], x_n);
                        }
                        else if (al.second == (*L)[x_n].back().vertex && Dis[x_n] > instance_graph.graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L,x, x_n)(x, x_n) + d_x) {
                            Dis[x_n] = instance_graph.graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L,x, x_n)(x, x_n) + d_x;
                            Q.emplace(Dis[x_n], x_n);
                        }
                    }
                }

        
                int h_c = al.second;
                if (Dis[x] == instance_graph.graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L,x, x_n)(x, al.first) + al.dis) 
				{
                    h_c = instance_graph.graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L,x, x_n)(al.first);
                }
                mtx_5952[x].lock();
                PPR_insert(*PPR, x, h_c, al.second);
                mtx_5952[x].unlock();
                mtx_5952[al.second].lock();
                PPR_insert(*PPR, al.second, h_c, x);
                mtx_5952[al.second].unlock();

                return 1;
            }
        }));
    }

    for (auto&& result : results_dynamic) {
        result.get();
    }
    std::vector<std::future<int>>().swap(results_dynamic);

}







void WeightDecreaseMaintenance_improv(graph_v_of_v_idealID& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	std::vector<affected_label> CL;
	WeightDecreaseMaintenance_improv_step1(v1, v2, w_new, &mm.L, &mm.PPR, &CL, pool_dynamic, results_dynamic);

	DIFFUSE(instance_graph, &mm.L, &mm.PPR, CL, pool_dynamic, results_dynamic);
}
