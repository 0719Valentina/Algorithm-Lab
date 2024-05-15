#pragma once

#include <build_in_progress/HL/dynamic/PLL_dynamic.h>
#include <map>

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
				results_dynamic.emplace_back(pool_dynamic.enqueue([it, v2, L, PPR, w_new, CL]
				{
					
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
	boost::heap::fibonacci_heap<PLL_dynamic_node_for_sp> Q;
	for (auto& cl : CL)
	{
		int u = cl.first;
		int v = cl.second;
		double du = cl.dis;
		
		std::vector<double>Dis(instance_graph.size(), -1.0);
		Dis[u] = du;
		/*u->du*/
		PLL_dynamic_node_for_sp node;
		node.vertex = u;
		node.priority_value = du;

		Q.push(node);

		//处理
		while (!Q.empty())
		{
			
			node=Q.top();
			int x = node.vertex;
			double dx = node.priority_value;
			Q.pop();

			vector<two_hop_label_v1>& L_x = (*L)[x];
			insert_sorted_two_hop_label(L_x, v, dx);
			for (const auto& neighbor : instance_graph[x])
			{
				int xn = neighbor.first;
				double w = neighbor.second;
				if (v < xn)
				{
					if (abs(Dis[xn]+1.0)<1e-5){
					Dis[xn] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(*L, v, xn);}
					double dnew = dx + w;

					if (Dis[xn] > dnew)
					{
						Dis[xn] = dnew;
						/*更新*/
						bool check=false;
						for(boost::heap::fibonacci_heap<PLL_dynamic_node_for_sp>::iterator it = Q.begin(); it != Q.end(); ++it){
							 if (it->vertex == xn) {
            					// 更新节点的优先级值
								boost::heap::fibonacci_heap<PLL_dynamic_node_for_sp>::handle_type handle = 
								Q.s_handle_from_iterator(it);
								(*handle).priority_value=Dis[xn];
            					// 调整 Fibonacci 堆以维持堆的性质
            					Q.update(handle);
								check=true;
            					break;
        					}
						}
						if(!check){
							PLL_dynamic_node_for_sp temp;
							temp.vertex = xn;
							temp.priority_value = Dis[xn];
							Q.push(temp);
						}
					}
					else
					{
						double Qxn=MAX_VALUE;
						bool check=false;
						boost::heap::fibonacci_heap<PLL_dynamic_node_for_sp>::iterator it = Q.begin();
						for(; it != Q.end(); ++it){
							 if (it->vertex == xn) {
            					Qxn=it->priority_value;
								check=true;
            					break;
        					}
						}
						auto result = search_sorted_two_hop_label((*L)[xn], v);
						int min = (Qxn > result) ? result : Qxn;
						if (result != MAX_VALUE && min > dnew){
							if(check){
								boost::heap::fibonacci_heap<PLL_dynamic_node_for_sp>::handle_type handle = 
								Q.s_handle_from_iterator(it);
								(*handle).priority_value=dnew;
            					// 调整 Fibonacci 堆以维持堆的性质
            					Q.update(handle);
							}
							else{
								PLL_dynamic_node_for_sp temp;
								temp.vertex = xn;
								temp.priority_value = dnew;
								Q.push(temp);
							}
						}
						auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, v, xn);
						if (query_result.second != v) {
							mtx_5952[xn].lock();
							PPR_insert(*PPR, xn, query_result.second, v);
							mtx_5952[xn].unlock();
						}
						if (query_result.second != xn) {
							mtx_5952[v].lock();
							PPR_insert(*PPR, v, query_result.second, xn);
							mtx_5952[v].unlock();
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
