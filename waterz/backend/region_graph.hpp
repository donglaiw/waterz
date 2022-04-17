#pragma once

#include "types.hpp"

#include <cstddef>
#include <iostream>
#include <map>

/**
 * Extract the region graph from a segmentation. Edges are annotated with the 
 * maximum affinity between the regions.
 *
 * @param aff [in]
 *              The affinity graph to read the affinities from.
 * @param seg [in]
 *              The segmentation.
 * @param max_segid [in]
 *              The highest ID in the segmentation.
 * @param start_zid [in]
 *              The starting z index. Default: 0
 * @param statisticsProvider [in]
 *              A statistics provider to update on-the-fly.
 * @param region_graph [out]
 *              A reference to a region graph to store the result.
 */
template<typename AG, typename V, typename StatisticsProviderType>
inline
void
get_region_graph(
		const AG& aff,
		const V& seg,
		std::size_t max_segid,
		std::size_t start_zid,
		StatisticsProviderType& statisticsProvider,
		RegionGraph<typename V::element>& rg) {

	typedef typename AG::element F;
	typedef typename V::element ID;
	typedef RegionGraph<ID> RegionGraphType;
	typedef typename RegionGraphType::EdgeIdType EdgeIdType;

	std::ptrdiff_t zdim = aff.shape()[1];
	std::ptrdiff_t ydim = aff.shape()[2];
	std::ptrdiff_t xdim = aff.shape()[3];

	// list of affinities between pairs of regions
	std::vector<std::map<ID, std::vector<F>>> affinities(max_segid+1);

	EdgeIdType e;
	std::size_t p[3];
    int cc = 0;
    // std::cout<<start_zid<<","<<zdim<<","<< ydim<<","<< xdim<<std::endl;
	for (p[0] = start_zid; p[0] < zdim; ++p[0])
		for (p[1] = 0; p[1] < ydim; ++p[1])
			for (p[2] = 0; p[2] < xdim; ++p[2]) {

				ID id1 = seg[p[0]][p[1]][p[2]];
                if (id1 == 0)
                    continue;
				statisticsProvider.addVoxel(id1, p[2], p[1], p[0]);

				for (int d = 0; d < 3; d++) {
					if (p[d] == 0)
						continue;

					ID id2 = seg[p[0]-(d==0)][p[1]-(d==1)][p[2]-(d==2)];
                    if (id2 == 0)
                        continue;

					if (id1 != id2) {
						auto mm = std::minmax(id1, id2);
						affinities[mm.first][mm.second].push_back(aff[d][p[0]][p[1]][p[2]]);
                        /*
                        if(mm.first==1 && mm.second==2){
                            std::cout<<"add: "<<mm.second<<"-"<<(int)aff[d][p[0]][p[1]][p[2]]<<","<<std::endl;
                        }
                        if (cc<30){
                            std::cout<<mm.first<<","<<mm.second<<","<< +aff[d][p[0]][p[1]][p[2]]<<std::endl;
                            cc+=1;
                        }
                        */
					}
				}
			}
	for (ID id1 = 1; id1 <= max_segid; ++id1) {
		for (const auto& p: affinities[id1]) {

			// p.first is ID
			// p.second is list of affiliated edges
			EdgeIdType e = rg.addEdge(id1, p.first);
			statisticsProvider.notifyNewEdge(e);

            // if (e == 0){std::cout<< (int)p.first<<"p"<<std::endl;}
			for (F affinity : p.second){
				statisticsProvider.addAffinity(e, affinity);
                // if ((int)id1 == 25356 && (int)p.first == 25604){std::cout<< e<<" "<<(int)affinity<<std::endl;}
                // if (e == 0){std::cout<< (int)affinity<<",";}
            }
        }
    }

	std::cout << "Region graph number of edges: " << rg.edges().size() << std::endl;
}


template<typename AG, typename V, typename StatisticsProviderType>
inline
void
get_region_graph_border(
		const AG& aff,
		const V& seg,
		std::size_t max_segid,
		StatisticsProviderType& statisticsProvider,
		RegionGraph<typename V::element>& rg) {

	typedef typename AG::element F;
	typedef typename V::element ID;
	typedef RegionGraph<ID> RegionGraphType;
	typedef typename RegionGraphType::EdgeIdType EdgeIdType;

	std::ptrdiff_t zdim = aff.shape()[1];
	std::ptrdiff_t ydim = aff.shape()[2];
	std::ptrdiff_t xdim = aff.shape()[3];

	// list of affinities between pairs of regions
	std::vector<std::map<ID, std::vector<F>>> affinities(max_segid+1);

	EdgeIdType e;
	std::size_t p[3];
    p[0] = 1;
		for (p[1] = 0; p[1] < ydim; ++p[1]){
			for (p[2] = 0; p[2] < xdim; ++p[2]) {

				ID id1 = seg[p[0]][p[1]][p[2]];
                if (id1 == 0)
                    continue;
				statisticsProvider.addVoxel(id1, p[2], p[1], p[0]);

                ID id2 = seg[p[0]-1][p[1]][p[2]];
                if (id2 == 0)
                    continue;
				statisticsProvider.addVoxel(id2, p[2], p[1], p[0]-1);

                if (id1 != id2) {
                    auto mm = std::minmax(id1, id2);
                    affinities[mm.first][mm.second].push_back(aff[0][p[0]][p[1]][p[2]]);
                }
			}
		}

	for (ID id1 = 1; id1 <= max_segid; ++id1) {
		for (const auto& p: affinities[id1]) {

			// p.first is ID
			// p.second is list of affiliated edges
			EdgeIdType e = rg.addEdge(id1, p.first);
			statisticsProvider.notifyNewEdge(e);

			for (F affinity : p.second)
				statisticsProvider.addAffinity(e, affinity);
        }
    }

	std::cout << "Region graph number of edges: " << rg.edges().size() << std::endl;
}

template<typename AG, typename V, typename StatisticsProviderType>
inline
void
get_region_graph_from_array(
        const V   num_edge,
		const AG& rg_score,
		const V&  rg_id1,
		const V&  rg_id2,
		StatisticsProviderType& statisticsProvider,
		RegionGraph<typename V::element>& rg) {

	typedef typename AG::element F;
	typedef typename V::element ID;
	typedef RegionGraph<ID> RegionGraphType;
	typedef typename RegionGraphType::EdgeIdType EdgeIdType;

	EdgeIdType e;

	for (ID i = 0; i < num_edge; ++i) {
			// p.first is ID
			// p.second is list of affiliated edges
			EdgeIdType e = rg.addEdge(rg_id1[i], rg_id2[i]);
			statisticsProvider.notifyNewEdge(e);
    }

	std::cout << "Region graph number of edges: " << rg.edges().size() << std::endl;
}
