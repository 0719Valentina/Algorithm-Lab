#pragma once

#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

void WeightDecreaseMaintenance_improv_step1(int v1, int v2, weightTYPE w_new, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR, std::vector<affected_label>* CL,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) 
{

	for (int sl = 0; sl < 2; sl++) {
		if (sl == 1) {
			swap(v1, v2);
		}
		//绗竴娆￠亶鍘哃(a) 绗簩娆￠亶鍘哃(b)
		for (auto it : (*L)[v1]) {
			//鍗硆(v)>=r(b) 浠庡皬鍒板ぇ鎺掑簭
			if (it.vertex <= v2) {
				results_dynamic.emplace_back(pool_dynamic.enqueue([it, v2, L, PPR, w_new, CL]
				{
					
					//璁＄畻Query璺濈
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
	for (auto& cl : CL)
	{
		int u = cl.first;
		int v = cl.second;
		double du = cl.dis;
		std::vector<double>Dis(instance_graph.size(), -1.0);
		Dis[u] = du;
		boost::heap::fibonacci_heap<PLL_dynamic_node_for_sp> Q;
		PLL_dynamic_node_for_sp node;
		node.vertex = u;
		node.priority_value = du;
		Q.push(node);

		while (!Q.empty())
		{

			node = Q.top();
			Q.pop();
			int x = node.vertex;
			weightTYPE dx = node.priority_value;
			vector<two_hop_label_v1>& L_x = (*L)[x];
			insert_sorted_two_hop_label(L_x, v, dx);
			for (const auto& neighbor : instance_graph[x])
			{
				int xn = neighbor.first;
				double w = neighbor.second;
				if (v < xn)
				{
					if (Dis[xn] == -1.0) Dis[xn] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(*L, xn, v);

					double dnew = dx + w;
					double Qxn = MAX_VALUE;
					PLL_dynamic_node_for_sp* check=nullptr;
					bool hasFind = false;
					for (auto n : Q)
					{
						if (n.vertex == xn)
						{
								check = &n;
								hasFind = true;
								break;
						}
					}

						if(Dis[xn] > dnew)
						{
							Dis[xn] = dnew;
							if (hasFind)
							{
								check->priority_value = Dis[xn];
							}else
							{
								PLL_dynamic_node_for_sp new_node;
								new_node.vertex = xn;
								new_node.priority_value = Dis[xn];
								Q.push(new_node);
							}
						}else
						{
							auto result = search_sorted_two_hop_label2((*L)[xn], v);
							if(hasFind)
							{
								Qxn=check->priority_value;
							}
							int min = (Qxn > result.first) ? result.first : Qxn;

							if (result.second != MAX_VALUE&& min > dnew)
							{
								if(hasFind)
								{
									check->priority_value = dnew;
								}
								else
								{
									PLL_dynamic_node_for_sp new_node;
									new_node.vertex = xn;
									new_node.priority_value = dnew;
									Q.push(new_node);
								}
							}
							
							auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, v, xn);
							if (query_result.second != x) {
								mtx_5952[xn].lock();
								PPR_insert(*PPR, xn, query_result.second, x);
								mtx_5952[xn].unlock();
							}
							if (query_result.second != xn) {
								mtx_5952[x].lock();
								PPR_insert(*PPR, x, query_result.second, xn);
								mtx_5952[x].unlock();
							}

						}
				}

				
			}

		}

	}
	for (auto&& result : results_dynamic){result.get();}
	std::vector<std::future<int>>().swap(results_dynamic);
}






void WeightDecreaseMaintenance_improv(graph_v_of_v_idealID& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic)
{
	std::vector<affected_label> CL;
	WeightDecreaseMaintenance_improv_step1(v1, v2, w_new, &mm.L, &mm.PPR, &CL, pool_dynamic, results_dynamic);
	DIFFUSE(instance_graph, &mm.L, &mm.PPR, CL, pool_dynamic, results_dynamic);
}
