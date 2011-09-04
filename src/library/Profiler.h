/* 
 * Copyright (C) 2008-2010  Daniel F. Schwarz
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PROFILER_H_
#define PROFILER_H_


#include <ctime>
#include <cassert>

class Profiler {
	public:
		/** 
		 * \brief Constructor
		 */
		Profiler ();

		/** 
		 * \brief Destructor
		 */
		~Profiler ();

		/** 
		 * \brief Start the stop watch
		 */
		void start();

		/** 
		 * \brief Stop the watch and add time difference to total time
		 */
		void stop();

		//! When did the profiler was starting the stop watch
		clock_t timeStamp;
		
		//! What is the total time that was consumed
		clock_t totalTime;

		//! Tick indicator
		bool isTicking;
};

#endif /*PROFILER_H_*/

