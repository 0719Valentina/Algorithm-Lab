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
	for(auto&cl:CL)
	{
		int u=cl.first;
		int v=cl.second;
		weightTYPE du=cl.dis;

		//初始化
		std::vector<weightTYPE>Dis(instance_graph.size(),-1);
		Dis[u]=du;
		std::priority_queue<std::pair<weightTYPE,int>>Q;
		Q.emplace(du,u);

		//处理
		while(!Q.empty())
		{
			weightTYPE dx=Q.top().first;
			int x=Q.top().second;
			Q.pop();
			(*L)[x][v].distance=dx;

			for(const auto&neighbor:instance_graph[x])
			{
				int xn=neighbor.first;
				weightTYPE w=neighbor.second;
				if(x>xn)
				{
					if(Dis[xn]==-1) Dis[xn]=graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(*L,xn,x);
					int dnew=dx+w;
					if(Dis[xn]>dnew)
					{
						Dis[xn]=dnew;
						Q.emplace(Dis[xn],xn);
					}else 
					{
						auto result = search_sorted_two_hop_label2((*L)[xn], x);
						int Qxn=0;
						std::priority_queue<std::pair<weightTYPE, int>> temp_q(Q); // 创建 Q 的副本
    					while (!temp_q.empty())
						{
        					auto [priority, vertex] = temp_q.top();
        					temp_q.pop();
        					if (vertex == xn) Qxn=priority;
    					}
						int min=(Qxn>result.first)?result.first:Qxn;
        				if (result.second != -1&&min>dnew)	Q.emplace(dnew,xn);

						 for (const auto& label_xn : L[xn])
						{
							for (const auto& label_x : L[x])
							{
								if (label_xn.vertex == label_x.vertex)
								{
									PPR_insert(*PPR,xn,label_x.vertex, x);
									PPR_insert(*PPR,x,label_x.vertex,xn);
								}
							}
						}
					}
				}

			}
			
		}

	}
	for (auto&& result : results_dynamic){
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
