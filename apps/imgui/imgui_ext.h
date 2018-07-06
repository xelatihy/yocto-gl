// Extensions for ImGUI by Fabio Pellacini

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

#ifndef _IMGUI_EXT_H_
#define _IMGUI_EXT_H_

#include "imgui.h"

#include <string>
#include <vector>
#include <functional>

// ImGui extensions
namespace ImGui {

// Check if active widgets
inline bool GetWidgetsActiveExt() {
    auto io = &ImGui::GetIO();
    return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

// Input text
inline bool InputText(const char* label, std::string* str) {
    char buf[4096];
    auto num = 0;
    for (auto c : *str) buf[num++] = c;
    buf[num] = 0;
    auto edited = InputText(label, buf, sizeof(buf));
    if (edited) *str = buf;
    return edited;
}

// Start selectable tree node
inline bool SelectableTreeNode(
    const char* lbl, void** selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
    if (*selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    auto open = ImGui::TreeNodeEx(content, node_flags, "%s", lbl);
    if (ImGui::IsItemClicked()) *selection = content;
    return open;
}

// Selectable tree leaf node
inline void SelectableTreeLeaf(
    const char* lbl, void*& selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;
    if (selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    ImGui::TreeNodeEx(content, node_flags, "%s", lbl);
    if (ImGui::IsItemClicked()) selection = content;
}

// Combo widget.
inline bool Combo(const char* lbl, int* idx_, int nitems,
    const std::function<const char*(int)>& label) {
    auto& idx = *idx_;
    if (!ImGui::BeginCombo(lbl, label(idx))) return false;
    auto old_idx = idx;
    for (auto i = 0; i < nitems; i++) {
        ImGui::PushID(i);
        if (ImGui::Selectable(label(i), idx == i)) idx = i;
        if (idx == i) ImGui::SetItemDefaultFocus();
        ImGui::PopID();
    }
    ImGui::EndCombo();
    return idx != old_idx;
}

// Combo widget.
inline bool Combo(const char* lbl, int* val_,
    const std::vector<std::string>& labels) {
    auto& val = *val_;
    if (!ImGui::BeginCombo(lbl, labels[val].c_str())) return false;
    auto old_val = val;
    for (auto i = 0; i < labels.size(); i++) {
        ImGui::PushID(i);
        if (ImGui::Selectable(labels[i].c_str(), val == i))
            val = i;
        if (val == i) ImGui::SetItemDefaultFocus();
        ImGui::PopID();
    }
    ImGui::EndCombo();
    return val != old_val;
}

// Combo widget.
inline bool Combo(const char* lbl, std::string* val_,
    const std::vector<std::string>& labels) {
    auto& val = *val_;
    if (!ImGui::BeginCombo(lbl, val.c_str())) return false;
    auto old_val = val;
    for (auto i = 0; i < labels.size(); i++) {
        ImGui::PushID(i);
        if (ImGui::Selectable(labels[i].c_str(), val == labels[i]))
            val = labels[i];
        if (val == labels[i]) ImGui::SetItemDefaultFocus();
        ImGui::PopID();
    }
    ImGui::EndCombo();
    return val != old_val;
}

// Combo widget
template <typename T>
inline bool Combo(
    const char* lbl, T** val_, const std::vector<T*>& vals, bool include_null) {
    auto& val = *val_;
    if (!ImGui::BeginCombo(lbl, (val) ? val->name.c_str() : "<none>"))
        return false;
    auto old_val = val;
    if (include_null) {
        ImGui::PushID(100000);
        if (ImGui::Selectable("<none>", val == nullptr)) val = nullptr;
        if (val == nullptr) ImGui::SetItemDefaultFocus();
        ImGui::PopID();
    }
    for (auto i = 0; i < vals.size(); i++) {
        ImGui::PushID(i);
        if (ImGui::Selectable(vals[i]->name.c_str(), val == vals[i]))
            val = vals[i];
        if (val == vals[i]) ImGui::SetItemDefaultFocus();
        ImGui::PopID();
    }
    ImGui::EndCombo();
    return val != old_val;
}

}  // namespace ImGui

#endif
