#pragma once

#include <queue>
#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

// ThreadPool
void SPREAD1(graph_v_of_v_idealID& instance_graph, vector<vector<two_hop_label_v1>>* L,
	std::vector<affected_label>& al1, std::vector<pair_label>* al2, ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	/*TO DO 2*/
	//遍历AL1
	for(auto it:al1){
		std::queue<std::pair<int,weightTYPE>> Q;
		Q.push(std::make_pair(it->first, it->dis));
		int v=it->second;
		while(!Q.empty()){
			std::pair<int,weightTYPE> temp=Q.front();
			Q.pop();
			int x=temp.first;
			weightTYPE dx=temp.second;
			L[x][v].distance=MAX_VALUE;
			al2.push_back(pair_label(x,v));
			//遍历x的邻接点
			int x_adj_size=ideal_graph_595[x].size();
			for(int i=0; i<x_adj_size; i++){
				int xn=ideal_graph_595[x][i].first;
				weightTYPE ec = ideal_graph_595[x][i].second;
				//r(v)>=r(xn)
				if(v<=xn){
					//if (𝑣, 𝑑𝑥 + 𝑤(𝑥, 𝑥𝑛 ) ) ∈ 𝐿(𝑥𝑛 ) then 𝑄𝑢𝑒𝑢𝑒.𝑝𝑢𝑠ℎ( (𝑥𝑛, 𝑑𝑥 + 𝑤(𝑥, 𝑥𝑛 ) ) )
					if(v<L[xn].size() && L[xn][v].distance==dx+ec)
					Q.push(std::make_pair(xn,dx+ec));
				}
			}
		}
	}
}

void SPREAD2(graph_v_of_v_idealID& instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR,
	std::vector<pair_label>& al2, std::vector<affected_label>* al3, ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	/*TO DO 3*/

}

void SPREAD3(graph_v_of_v_idealID& instance_graph, vector<vector<two_hop_label_v1>>* L, PPR_type* PPR, std::vector<affected_label>& al3,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	/*TO DO 4*/
	for(auto it:al3){
		int u=it->first;
		int v=it->second;
		weightTYPE du=it->dis;

		auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, u, v);

		//求query的值

			if(query_result.first <= du){
				if (query_result.second != v) {
						mtx_5952[v2].lock();
						PPR_insert(*PPR, u, query_result.second, v);
						mtx_5952[v2].unlock();
					}
				if (query_result.second != u) {
						mtx_5952[it.vertex].lock();
						PPR_insert(*PPR, v, query_result.second, u);
						mtx_5952[it.vertex].unlock();
					}
				continue;
			}

		//初始化dis数组 
		std::vector<int> DIS;
		// typedef std::vector<std::vector<std::pair<int, double>>> graph_v_of_v_idealID;
		int v_size=instance_graph.size();
		for(int i=0; i<v_size; i++){
			if(i==u)
			DIS[i]=du;
			else
			DIS[i]=-1;
		}

		//初始化Q
		std::queue<std::pair<int,weightTYPE>> Q;
		Q.push(std::make_pair(u, du));

		while(!Q.empty()){
			std::pair<int,weightTYPE> temp=Q.front();
			Q.pop();
			int x=temp.first;
			weightTYPE dx=temp.second;

			if(dx < L[x][v].distance)
			L[x][v].distance=dx;

			//遍历x的邻接点
			int x_adj_size=ideal_graph_595[x].size();
			for(int i=0; i<x_adj_size; i++){
				int xn=ideal_graph_595[x][i].first;
				weightTYPE ec = ideal_graph_595[x][i].second;
				//r(v)>=r(xn)
				if(v<=xn){
					if(DIS[xn]==-1)
					DIS[xn]=graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, xn, v);
					
					if(DIS[xn] > dx+ec){
						DIS[xn]=dx+ec;
						//update
					}
				}
			}
		}
	}
	}

//G、L、PPR、a、b、w0 v1->a v2->b mm.L->v+d 
void WeightIncreaseMaintenance_improv(graph_v_of_v_idealID& instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1& mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new,
	ThreadPool& pool_dynamic, std::vector<std::future<int>>& results_dynamic) {

	std::vector<affected_label> al1, al3;
	std::vector<pair_label> al2;

	/*it's slow to paralize the following part*/
	for (auto it : mm.L[v1]) {
		if (it.vertex <= v2 && abs(search_sorted_two_hop_label(mm.L[v2], it.vertex) - it.distance - w_old) < 1e-5) {
			al1.push_back(affected_label(v2, it.vertex, it.distance + w_old));
		}
	}
	for (auto it : mm.L[v2]) {
		if (it.vertex <= v1 && abs(search_sorted_two_hop_label(mm.L[v1], it.vertex) - it.distance - w_old) < 1e-5) {
			al1.push_back(affected_label(v1, it.vertex, it.distance + w_old));
		}
	}

	//cout << "al1.size() " << al1.size() << endl;

	SPREAD1(instance_graph, &mm.L, al1, &al2, pool_dynamic, results_dynamic);
	SPREAD2(instance_graph, &mm.L, &mm.PPR, al2, &al3, pool_dynamic, results_dynamic);
	SPREAD3(instance_graph, &mm.L, &mm.PPR, al3, pool_dynamic, results_dynamic);

	//for (auto it : al2) {
	//	cout << "al2 " << it.first << " " << it.second << endl;
	//}
	//for (auto it : al3) {
	//	cout << "al3 " << it.first << " " << it.second << " " << it.dis << endl;
	//}
}

