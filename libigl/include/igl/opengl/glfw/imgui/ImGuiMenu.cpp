// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Jérémie Dumas <jeremie.dumas@ens-lyon.org>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
////////////////////////////////////////////////////////////////////////////////
#include "ImGuiMenu.h"
#include <igl/project.h>
#include <imgui/imgui.h>
#include <imgui_impl_glfw_gl3.h>
#include <imgui_fonts_droid_sans.h>
#include <GLFW/glfw3.h>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace igl
{
namespace opengl
{
namespace glfw
{
namespace imgui
{

IGL_INLINE void ImGuiMenu::init(igl::opengl::glfw::Viewer *_viewer)
{
  ViewerPlugin::init(_viewer);
  // Setup ImGui binding
  if (_viewer)
  {
    if (context_ == nullptr)
    {
      context_ = ImGui::CreateContext();
    }
    ImGui_ImplGlfwGL3_Init(viewer->window, false);
    ImGui::GetIO().IniFilename = nullptr;
    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    style.FrameRounding = 5.0f;
    reload_font();
  }
}

IGL_INLINE void ImGuiMenu::reload_font(int font_size)
{
  hidpi_scaling_ = hidpi_scaling();
  pixel_ratio_ = pixel_ratio();
  ImGuiIO& io = ImGui::GetIO();
  io.Fonts->Clear();
  io.Fonts->AddFontFromMemoryCompressedTTF(droid_sans_compressed_data,
    droid_sans_compressed_size, font_size * hidpi_scaling_);
  io.FontGlobalScale = 1.0 / pixel_ratio_;
}

IGL_INLINE void ImGuiMenu::shutdown()
{
  // Cleanup
  ImGui_ImplGlfwGL3_Shutdown();
  ImGui::DestroyContext(context_);
  context_ = nullptr;
}

IGL_INLINE bool ImGuiMenu::pre_draw()
{
  glfwPollEvents();

  // Check whether window dpi has changed
  float scaling = hidpi_scaling();
  if (std::abs(scaling - hidpi_scaling_) > 1e-5)
  {
    reload_font();
    ImGui_ImplGlfwGL3_InvalidateDeviceObjects();
  }

  ImGui_ImplGlfwGL3_NewFrame();
  return false;
}

IGL_INLINE bool ImGuiMenu::post_draw()
{
  draw_menu();
  ImGui::Render();
  return false;
}

IGL_INLINE void ImGuiMenu::post_resize(int width, int height)
{
  if (context_)
  {
    ImGui::GetIO().DisplaySize.x = float(width);
    ImGui::GetIO().DisplaySize.y = float(height);
  }
}

// Mouse IO
IGL_INLINE bool ImGuiMenu::mouse_down(int button, int modifier)
{
  ImGui_ImplGlfwGL3_MouseButtonCallback(viewer->window, button, GLFW_PRESS, modifier);
  return ImGui::GetIO().WantCaptureMouse;
}

IGL_INLINE bool ImGuiMenu::mouse_up(int button, int modifier)
{
  return ImGui::GetIO().WantCaptureMouse;
}

IGL_INLINE bool ImGuiMenu::mouse_move(int mouse_x, int mouse_y)
{
  return ImGui::GetIO().WantCaptureMouse;
}

IGL_INLINE bool ImGuiMenu::mouse_scroll(float delta_y)
{
  ImGui_ImplGlfwGL3_ScrollCallback(viewer->window, 0.f, delta_y);
  return ImGui::GetIO().WantCaptureMouse;
}

// Keyboard IO
IGL_INLINE bool ImGuiMenu::key_pressed(unsigned int key, int modifiers)
{
  ImGui_ImplGlfwGL3_CharCallback(nullptr, key);
  return ImGui::GetIO().WantCaptureKeyboard;
}

IGL_INLINE bool ImGuiMenu::key_down(int key, int modifiers)
{
  ImGui_ImplGlfwGL3_KeyCallback(viewer->window, key, 0, GLFW_PRESS, modifiers);
  return ImGui::GetIO().WantCaptureKeyboard;
}

IGL_INLINE bool ImGuiMenu::key_up(int key, int modifiers)
{
  ImGui_ImplGlfwGL3_KeyCallback(viewer->window, key, 0, GLFW_RELEASE, modifiers);
  return ImGui::GetIO().WantCaptureKeyboard;
}

// Draw menu
IGL_INLINE void ImGuiMenu::draw_menu()
{
  // Text labels
  draw_labels_window();

  // Viewer settings
  float menu_width = 180.f * menu_scaling();
  ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
  ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
  ImGui::SetNextWindowSizeConstraints(ImVec2(menu_width, -1.0f), ImVec2(menu_width, -1.0f));
  bool _viewer_menu_visible = true;
  ImGui::Begin(
      "Menu", &_viewer_menu_visible,
      ImGuiWindowFlags_NoSavedSettings
      | ImGuiWindowFlags_AlwaysAutoResize
  );
  ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.4f);
  if (draw_viewer_menu_func) { draw_viewer_menu_func(); }
  else { draw_viewer_menu(); }
  ImGui::PopItemWidth();
  ImGui::End();

  // Other windows
  if (draw_custom_window_func) { draw_custom_window_func(); }
  else { draw_custom_window(); }
}

IGL_INLINE void ImGuiMenu::draw_viewer_menu()
{
 
}

IGL_INLINE void ImGuiMenu::draw_labels_window()
{
  // Text labels
  ImGui::SetNextWindowPos(ImVec2(0,0), ImGuiSetCond_Always);
  ImGui::SetNextWindowSize(ImGui::GetIO().DisplaySize, ImGuiSetCond_Always);
  bool visible = true;
  ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0,0,0,0));
  ImGui::Begin("ViewerLabels", &visible,
      ImGuiWindowFlags_NoTitleBar
      | ImGuiWindowFlags_NoResize
      | ImGuiWindowFlags_NoMove
      | ImGuiWindowFlags_NoScrollbar
      | ImGuiWindowFlags_NoScrollWithMouse
      | ImGuiWindowFlags_NoCollapse
      | ImGuiWindowFlags_NoSavedSettings
      | ImGuiWindowFlags_NoInputs);
  for (const auto & data : viewer->data_list)
  {
    draw_labels(data);
  }
  ImGui::End();
  ImGui::PopStyleColor();
}

IGL_INLINE void ImGuiMenu::draw_labels(const igl::opengl::ViewerData &data)
{
  if (data.show_vertid)
  {
    for (int i = 0; i < data.V.rows(); ++i)
    {
      draw_text(data.V.row(i), data.V_normals.row(i), std::to_string(i));
    }
  }

  if (data.show_faceid)
  {
    for (int i = 0; i < data.F.rows(); ++i)
    {
      Eigen::RowVector3d p = Eigen::RowVector3d::Zero();
      for (int j = 0; j < data.F.cols(); ++j)
      {
        p += data.V.row(data.F(i,j));
      }
      p /= (double) data.F.cols();

      draw_text(p, data.F_normals.row(i), std::to_string(i));
    }
  }

  if (data.labels_positions.rows() > 0)
  {
    for (int i = 0; i < data.labels_positions.rows(); ++i)
    {
      draw_text(data.labels_positions.row(i), Eigen::Vector3d(0.0,0.0,0.0),
        data.labels_strings[i]);
    }
  }
}

IGL_INLINE void ImGuiMenu::draw_text(Eigen::Vector3d pos, Eigen::Vector3d normal, const std::string &text)
{
  Eigen::Matrix4f view_matrix = viewer->core.view * viewer->core.model;
  pos += normal * 0.005f * viewer->core.object_scale;
  Eigen::Vector3f coord = igl::project(Eigen::Vector3f(pos.cast<float>()),
    view_matrix, viewer->core.proj, viewer->core.viewport);

  // Draw text labels slightly bigger than normal text
  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->AddText(ImGui::GetFont(), ImGui::GetFontSize() * 1.2,
      ImVec2(coord[0]/pixel_ratio_, (viewer->core.viewport[3] - coord[1])/pixel_ratio_),
      ImGui::GetColorU32(ImVec4(0, 0, 10, 255)),
      &text[0], &text[0] + text.size());
}

IGL_INLINE float ImGuiMenu::pixel_ratio()
{
  // Computes pixel ratio for hidpi devices
  int buf_size[2];
  int win_size[2];
  GLFWwindow* window = glfwGetCurrentContext();
  glfwGetFramebufferSize(window, &buf_size[0], &buf_size[1]);
  glfwGetWindowSize(window, &win_size[0], &win_size[1]);
  return (float) buf_size[0] / (float) win_size[0];
}

IGL_INLINE float ImGuiMenu::hidpi_scaling()
{
  // Computes scaling factor for hidpi devices
  float xscale, yscale;
  GLFWwindow* window = glfwGetCurrentContext();
  glfwGetWindowContentScale(window, &xscale, &yscale);
  return 0.5 * (xscale + yscale);
}

} // end namespace
} // end namespace
} // end namespace
} // end namespace