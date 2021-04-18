#pragma once

#include "query/queries.hpp"
#include <vector>

namespace pisa {

template <typename Index>
[[nodiscard]] auto make_cursors(Index const& index, Query query)
{
    auto terms = query.terms;
    remove_duplicate_terms(terms);
    using cursor = typename Index::document_enumerator;

    std::vector<cursor> cursors;
    cursors.reserve(terms.size());
    std::transform(terms.begin(), terms.end(), std::back_inserter(cursors), [&](auto&& term) {
        /*double tick = get_time_usecs();
        cursor indexterm = index[term];
        double elapsed_secs = (get_time_usecs() - tick) / 1000;
        std::cout<<"cursor - "<<elapsed_secs<<" miliseconds"<<std::endl;
        return indexterm;*/
        return index[term];
    });

    return cursors;
}

}  // namespace pisa
