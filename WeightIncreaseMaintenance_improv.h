#pragma once

#include <queue>
#include <build_in_progress/HL/dynamic/PLL_dynamic.h>

void SPREAD1(graph_v_of_v_idealID &instance_graph, vector<vector<two_hop_label_v1>> *L,
			 std::vector<affected_label> &al1, std::vector<pair_label> *al2, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{

	/*TO DO 2*/
	for (auto it : al1)
	{
		std::queue<std::pair<int, weightTYPE>> Q;
		Q.push(std::make_pair(it->first, it->dis));
		int v = it->second;
		while (!Q.empty())
		{
			std::pair<int, weightTYPE> temp = Q.front();
			Q.pop();
			int x = temp.first;
			weightTYPE dx = temp.second;
			L[x][v].distance = MAX_VALUE;
			al2.push_back(pair_label(x, v));
			// 遍历x的邻接点
			//int x_adj_size = ideal_graph_595[x].size();
			for (const auto &neighbor : instance_graph[x]) // for (int i = 0; i < x_adj_size; i++)
			{
				int xn = neighbor.first;
				weightTYPE ec = neighbor.second;
				// r(v)>=r(xn)
				if (v <= xn)
				{
					// if (𝑣, 𝑑𝑥 + 𝑤(𝑥, 𝑥𝑛 ) ) ∈ 𝐿(𝑥𝑛 ) then 𝑄𝑢𝑒𝑢𝑒.𝑝𝑢𝑠ℎ( (𝑥𝑛, 𝑑𝑥 + 𝑤(𝑥, 𝑥𝑛 ) ) )
					auto search_result = search_sorted_two_hop_label((*L)[xn], v);
					if (search_result == dx + ec)
						Q.push(std::make_pair(xn, dx + ec));
				}
			}
		}
	}
}

void SPREAD2(graph_v_of_v_idealID &instance_graph, vector<vector<two_hop_label_v1>> *L, PPR_type *PPR,
			 std::vector<pair_label> &al2, std::vector<affected_label> *al3, ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{

	for (auto it : al2)
	{
		int x = it.first;
		int y = it.second;
		// weightTYPE w = it.dis;
		for (auto t : (PPR[x][y] || y)) // If 𝑡 ∈ 𝑃𝑃𝑅[𝑥, 𝑦] ∪ 𝑦
		{
			// if 𝑟 (𝑡) > 𝑟 (𝑥 )
			if (t < x)
			{
				// 在xn中循环找到最小值
				weightTYPE d1x_t = MAX_VALUE;				   // 初始化无穷大
				for (const auto &neighbor : instance_graph[x]) // for (int i = 0; i < x_adj_size; i++)
				{
					int xn = neighbor.first;
					auto search_result = search_sorted_two_hop_label((*L)[xn], t);
					weightTYPE ec = neighbor.second;
					if (d1x_t > ec + search_result)
					{
						d1x_t = ec + search_result;
					}
				}

				// if 𝑄𝑢𝑒𝑟𝑦(𝑥, 𝑡, 𝐿) > 𝑑1(𝑥, 𝑡) then 𝐴𝐿3.𝑝𝑢𝑠ℎ( (𝑥, 𝑡, 𝑑1(𝑥, 𝑡) ) )
				if (graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(*L, x, t) > d1x_t)
				{
					al3.push_back(affected_label(x, t, d1x_t));
				}
				else // else 𝑃𝑃𝑅[𝑥, ℎ𝑐 ].𝑝𝑢𝑠ℎ(𝑡), 𝑃𝑃𝑅[𝑡, ℎ𝑐 ].𝑝𝑢𝑠ℎ(𝑥 )
				{
					auto query_result2 = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, x, t);

					PPR_insert(*PPR, x, query_result2.second, t);
					PPR_insert(*PPR, t, query_result2.second, x);

					/*
					auto query_result2 = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, x, t);
					if (query_result2.second != t)
					{
						mtx_5952[x].lock();
						PPR_insert(*PPR, x, query_result2.second, t);
						mtx_5952[x].unlock();
					}
					if (query_result2.second != x)
					{
						mtx_5952[t].lock();
						PPR_insert(*PPR, t, query_result2.second, x);
						mtx_5952[t].unlock();
					}
					*/
				}
			}
			else if (t > x)
			{
				// 在xn中循环找到最小值
				weightTYPE d1t_x = MAX_VALUE; // 初始化无穷大
				// int t_adj_size = ideal_graph_595[t].size();
				for (const auto &neighbor : instance_graph[x]) // for (int i = 0; i < x_adj_size; i++)
				{
					int tn = neighbor.first;
					auto search_result = search_sorted_two_hop_label((*L)[tn], x);
					weightTYPE ec = neighbor.second;
					if (d1t_x > ec + search_result)
					{
						d1t_x = ec + search_result;
					}
				}

				// if 𝑄𝑢𝑒𝑟𝑦(𝑥, 𝑡, 𝐿) > 𝑑1(𝑥, 𝑡) then 𝐴𝐿3.𝑝𝑢𝑠ℎ( (𝑥, 𝑡, 𝑑1(𝑥, 𝑡) ) )
				if (graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(*L, t, x) > d1t_x)
				{
					al3.push_back(affected_label(t, x, d1t_x));
				}
				else // else 𝑃𝑃𝑅[𝑥, ℎ𝑐 ].𝑝𝑢𝑠ℎ(𝑡), 𝑃𝑃𝑅[𝑡, ℎ𝑐 ].𝑝𝑢𝑠ℎ(𝑥 )
				{
					auto query_result2 = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, t, x);

					PPR_insert(*PPR, t, query_result2.second, x);
					PPR_insert(*PPR, x, query_result2.second, t);

					/*
					auto query_result2 = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, t, x);
					if (query_result2.second != x)
					{
						mtx_5952[t].lock();
						PPR_insert(*PPR, t, query_result2.second, x);
						mtx_5952[t].unlock();
					}
					if (query_result2.second != t)
					{
						mtx_5952[x].lock();
						PPR_insert(*PPR, x, query_result2.second, t);
						mtx_5952[x].unlock();
					}
					*/
				}
			}
		}
	}
}

void SPREAD3(graph_v_of_v_idealID &instance_graph, vector<vector<two_hop_label_v1>> *L, PPR_type *PPR, std::vector<affected_label> &al3,
			 ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{

	/*TO DO 4*/
	for (auto it : al3)
	{
		int u = it->first;
		int v = it->second;
		weightTYPE du = it->dis;

		auto query_result = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, u, v);

		// 求query的值

		if (query_result.first <= du)
		{
			if (query_result.second != v)
			{
				mtx_5952[u].lock();
				PPR_insert(*PPR, u, query_result.second, v);
				mtx_5952[u].unlock();
			}
			if (query_result.second != u)
			{
				mtx_5952[v].lock();
				PPR_insert(*PPR, v, query_result.second, u);
				mtx_5952[v].unlock();
			}
			continue;
		}

		// 初始化dis数组
		std::vector<int> DIS;
		int v_size = instance_graph.size();
		for (int i = 0; i < v_size; i++)
		{
			if (i == u)
				DIS[i] = du;
			else
				DIS[i] = -1;
		}

		// 初始化Q 斐波那契堆 最小堆
		boost::heap::fibonacci_heap<PLL_dynamic_node_for_sp> Q;
		PLL_dynamic_node_for_sp node;
		node.vertex = u;
		node.priority_value = du;

		Q.push(node);

		while (Q.size() > 0)
		{
			node = Q.top();
			Q.pop();

			int x = node.vertex;
			weightTYPE dx = node.priority_value;

			// 取Min值
			if (dx < L[x][v].distance)
				L[x][v].distance = dx;

			// 遍历x的邻接点
			// int x_adj_size = ideal_graph_595[x].size();
			for (const auto &neighbor : instance_graph[x]) // for (int i = 0; i < x_adj_size; i++)
			{
				int xn = neighbor.first;
				weightTYPE ec = neighbor.second;
				// r(v)>=r(xn)
				if (v <= xn)
				{
					if (DIS[xn] == -1)
						DIS[xn] = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(*L, xn, v); // query_result is {distance, common hub};

					if (DIS[xn] > dx + ec)
					{
						DIS[xn] = dx + ec;
						// update 查看Q中是否有xn
						auto it = std::find_if(Q.begin(), Q.end(), [xn](const PLL_dynamic_node_for_sp &node)
											   { return node.vertex == xn; });
						if (it != Q.end())
						{
							it->priority_value = DIS[xn];
							Q.update(it);
						}
						else
						{
							PLL_dynamic_node_for_sp new_node;
							new_node.vertex = xn;
							new_node.priority_value = DIS[xn];
							Q.push(new_node);
						}
					}
					else
					{
						auto query_result2 = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc2(*L, xn, v);
						if (query_result2.second != v)
						{
							mtx_5952[xn].lock();
							PPR_insert(*PPR, xn, query_result2.second, v);
							mtx_5952[xn].unlock();
						}
						if (query_result2.second != xn)
						{
							mtx_5952[v].lock();
							PPR_insert(*PPR, v, query_result2.second, xn);
							mtx_5952[v].unlock();
						}
					}
				}
			}
		}
	}
}

void WeightIncreaseMaintenance_improv(graph_v_of_v_idealID &instance_graph, graph_hash_of_mixed_weighted_two_hop_case_info_v1 &mm, int v1, int v2, weightTYPE w_old, weightTYPE w_new,
									  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic)
{

	std::vector<affected_label> al1, al3;
	std::vector<pair_label> al2;

	/*it's slow to paralize the following part*/
	for (auto it : mm.L[v1])
	{
		if (it.vertex <= v2 && abs(search_sorted_two_hop_label(mm.L[v2], it.vertex) - it.distance - w_old) < 1e-5)
		{
			al1.push_back(affected_label(v2, it.vertex, it.distance + w_old));
		}
	}
	for (auto it : mm.L[v2])
	{
		if (it.vertex <= v1 && abs(search_sorted_two_hop_label(mm.L[v1], it.vertex) - it.distance - w_old) < 1e-5)
		{
			al1.push_back(affected_label(v1, it.vertex, it.distance + w_old));
		}
	}

	// cout << "al1.size() " << al1.size() << endl;

	SPREAD1(instance_graph, &mm.L, al1, &al2, pool_dynamic, results_dynamic);
	SPREAD2(instance_graph, &mm.L, &mm.PPR, al2, &al3, pool_dynamic, results_dynamic);
	SPREAD3(instance_graph, &mm.L, &mm.PPR, al3, pool_dynamic, results_dynamic);

	// for (auto it : al2) {
	//	cout << "al2 " << it.first << " " << it.second << endl;
	// }
	// for (auto it : al3) {
	//	cout << "al3 " << it.first << " " << it.second << " " << it.dis << endl;
	// }
}
