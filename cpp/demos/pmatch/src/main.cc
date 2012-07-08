/*
 * LibSRVF - a shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2012  Daniel Robinson
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */
#include "partialmatch.h"
#include "paretoset.h"
#include "matchview.h"
#include "ui.h"

#include <matrix.h>
#include <plf.h>
#include <srvf.h>
#include <qmap.h>
#include <fileio.h>

#include <FL/Fl.H>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

#define DEFAULT_GRID_WIDTH 31
#define DEFAULT_GRID_HEIGHT 31


static void do_usage(const char *progname)
{
  std::cout << "USAGE: " << progname << " [OPTIONS] file1 file2";
  std::cout << "where" << std::endl;
  std::cout << "  file1 contains samples for the first curve" << std::endl;
  std::cout << "  file2 contains samples for the second curve" << std::endl;
  std::cout << "Recognized options are:" << std::endl;
  std::cout << "  -W n\t\tuse grid width n" << std::endl;
  std::cout << "  -H n\t\tuse grid height n" << std::endl;
  std::cout << "  -r\t\toptimize over rotations" << std::endl;
  std::cout << "  -o outfile\t\twrite matches to outfile" << std::endl;
  std::cout << "  -h\t\tshow this message" << std::endl;
}

int main( int argc, char **argv ){
  size_t grid_width = DEFAULT_GRID_WIDTH;
  size_t grid_height = DEFAULT_GRID_HEIGHT;
  bool do_rotations = false;
  const char *output_filename = "matches.mat";
  int opt;

  while( (opt=getopt(argc, argv, "hH:o:rW:")) != -1 ){
    switch( opt ){
      case 'h':
        do_usage(argv[0]);
        return 0;
      case 'H':
        grid_height = atoi(optarg);
        break;
      case 'o':
        output_filename = optarg;
        break;
      case 'r':
        do_rotations = true;
        break;
      case 'W':
        grid_width = atoi(optarg);
        break;
      case '?':
        do_usage(argv[0]);
        return -1;
    }
  }
  if (argc - optind != 2) {
    do_usage(argv[0]);
    return -1;
  }

  std::ifstream ifs1(argv[optind]);
  std::ifstream ifs2(argv[optind+1]);

  std::vector<srvf::Matrix> F1data = srvf::io::load_csv ( ifs1, ' ', '\n' );
  std::vector<srvf::Matrix> F2data = srvf::io::load_csv ( ifs2, ' ', '\n' );

  ifs1.close();
  ifs2.close();

  srvf::Pointset F1samps(F1data[0], srvf::Pointset::POINT_PER_COLUMN);
  srvf::Pointset F2samps(F2data[0], srvf::Pointset::POINT_PER_COLUMN);

  srvf::Plf F1(F1samps);
  srvf::Plf F2(F2samps);

  double L1 = F1.arc_length();
  double L2 = F2.arc_length();
  double L = std::min(L1, L2);
  F1.scale( 1.0 / L );
  F2.scale( 1.0 / L );

  srvf::Srvf Q1 = srvf::plf_to_srvf(F1);
  srvf::Srvf Q2 = srvf::plf_to_srvf(F2);

  srvf::pmatch::ParetoSet S = srvf::pmatch::find_matches (
    Q1, Q2, false, grid_width, grid_height);

  std::ofstream ofs(output_filename);
  S.save_csv(ofs);
  ofs.close();

  srvf::SuperimposedPlot plot;
  plot.insert(F1, srvf::Color(0.0, 0.0, 1.0));
  plot.insert(F2, srvf::Color(1.0, 0.0, 0.0));

  UserInterface ui;

  Fl_Window *win = ui.make_window();
  ui.match_view->set_plot(plot);
  ui.match_view->set_matches(S);
  ui.slider_match_length->range(0.0, (double)(S.nbuckets()-1));
  ui.slider_match_length->step(1.0);
  ui.slider_match_length->value(1.0);

  win->show();
  Fl::run();

  return 0;
}

