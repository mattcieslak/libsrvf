// generated by Fast Light User Interface Designer (fluid) version 1.0300

#include "ui.h"

void UserInterface::cb_slider_match_length_i(Fl_Slider* o, void*) {
  size_t idx=(size_t)round(o->value());
match_view->set_match(idx,0);
}
void UserInterface::cb_slider_match_length(Fl_Slider* o, void* v) {
  ((UserInterface*)(o->parent()->parent()->user_data()))->cb_slider_match_length_i(o,v);
}

Fl_Window* UserInterface::make_window() {
  Fl_Window* w;
  { Fl_Window* o = new Fl_Window(695, 555, "Partial Match");
    w = o;
    o->box(FL_FLAT_BOX);
    o->color(FL_BACKGROUND_COLOR);
    o->selection_color(FL_BACKGROUND_COLOR);
    o->labeltype(FL_NO_LABEL);
    o->labelfont(0);
    o->labelsize(14);
    o->labelcolor(FL_FOREGROUND_COLOR);
    o->user_data((void*)(this));
    o->align(Fl_Align(FL_ALIGN_TOP));
    o->when(FL_WHEN_RELEASE);
    { Fl_Pack* o = new Fl_Pack(0, 0, 695, 555);
      { match_view = new MatchView(5, 3, 685, 494);
        match_view->box(FL_NO_BOX);
        match_view->color(FL_BACKGROUND_COLOR);
        match_view->selection_color(FL_BACKGROUND_COLOR);
        match_view->labeltype(FL_NORMAL_LABEL);
        match_view->labelfont(0);
        match_view->labelsize(14);
        match_view->labelcolor(FL_FOREGROUND_COLOR);
        match_view->align(Fl_Align(FL_ALIGN_CENTER));
        match_view->when(FL_WHEN_RELEASE);
      } // MatchView* match_view
      { slider_match_length = new Fl_Slider(5, 505, 685, 30, "Match Length");
        slider_match_length->type(3);
        slider_match_length->box(FL_FLAT_BOX);
        slider_match_length->color((Fl_Color)55);
        slider_match_length->maximum(20);
        slider_match_length->step(1);
        slider_match_length->value(1);
        slider_match_length->slider_size(1);
        slider_match_length->callback((Fl_Callback*)cb_slider_match_length);
        slider_match_length->align(Fl_Align(290));
        Fl_Group::current()->resizable(slider_match_length);
      } // Fl_Slider* slider_match_length
      o->end();
    } // Fl_Pack* o
    o->end();
    o->resizable(o);
  } // Fl_Window* o
  return w;
}
