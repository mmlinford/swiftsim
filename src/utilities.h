/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef SWIFT_UTILITIES_H
#define SWIFT_UTILITIES_H

/**
 * @brief Search for a value in a monotonically increasing array to find the
 *      index such that array[index] < value < array[index + 1]
 *
 * @param x The value to find
 * @param array The array to search
 * @param n The length of the array
 *
 * Return -1 and n for x below and above the array edge values respectively.
 */
INLINE static int find_value_in_monotonic_array(
    const float x, const float *array, const int n) {

    int index_mid, index_low = 0, index_high = n;

    // Until array[index_low] < x < array[index_high=index_low+1]
    while (index_high - index_low > 1) {
        index_mid = (index_high + index_low) / 2.f;  // Middle index

        // Replace the low or high index with the middle
        if (array[index_mid] <= x)
            index_low = index_mid;
        else
            index_high = index_mid;
    }

    // Set index with the found index_low or an error value if outside the array
    if (x < array[0])
        return -1;
    else if (array[n-1] <= x)
        return n;
    else
        return index_low;
}

#endif /* SWIFT_UTILITIES_H */
