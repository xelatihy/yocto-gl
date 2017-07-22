///
/// YOCTO_DOCGEN: Extract documentation from yocto.
///
/// This code is realsed only as a demonstration for others.
///

//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#include "../yocto/yocto_utils.h"
using namespace yu::operators;
using namespace yu::string;
using namespace yu::file;

// clang-format off
auto filenames = std::vector<std::string>{
    "yocto/yocto_bvh.h",
    "yocto/yocto_gltf.h",
    "yocto/yocto_glu.h",
    "yocto/yocto_gui.h",
    "yocto/yocto_img.h",
    "yocto/yocto_math.h",
    "yocto/yocto_obj.h",
    "yocto/yocto_sym.h",
    "yocto/yocto_utils.h",
    "yocto/yocto_trace.h",
};
// clang-format on

struct item {
    std::string name;
    std::string decl;
    std::string comment;
    std::string type;
    std::vector<item> children;
};

std::string make_doc(const std::string& cpp) {
    // comment blocks
    auto items = std::vector<item>{};
    auto cur_item = (item*)nullptr;
    bool first = true, indented = false, enum_last = false;
    for (auto line : splitlines(cpp, true)) {
        if (cur_item) {
            if (contains(line, "///")) {
                cur_item->comment += line;
            } else if (first) {
                cur_item = nullptr;
                first = false;
            } else {
                cur_item->decl += line;
                if (contains(line, ";") || contains(line, "{")) {
                    cur_item = nullptr;
                }
                if (enum_last && contains(line, ",")) { cur_item = nullptr; }
            }
        } else {
            if (!contains(line, "///")) continue;
            if (startswith(line, "    ")) {
                items.back().children.push_back({});
                cur_item = &items.back().children.back();
                cur_item->comment += line;
                indented = true;
                enum_last = contains(items.back().decl, "enum ");
            } else {
                items.push_back({});
                cur_item = &items.back();
                cur_item->comment += line;
                indented = false;
                enum_last = false;
            }
        }
    }

    auto clean_comment = [](std::string comment) {
        if (startswith(comment, "    ")) {
            comment = replace(comment, "    ///", "");
            comment = replace(comment, "    /// ", "");
            comment = strip(comment);
        } else {
            comment = replace(comment, "/// ", "");
            comment = replace(comment, "///", "");
            comment = strip(comment);
        }
        return comment;
    };

    // fix type
    for (auto& item : items) {
        if (item.decl == "") {
            item.comment = clean_comment(item.comment);
        } else if (contains(item.decl, "namespace ")) {
            item.type = "Namespace";
            item.name = item.decl;
            item.name = replace(item.name, "namespace ", "");
            item.name = replace(item.name, "{", "");
            item.name = strip(item.name);
            item.decl = "";
            item.comment = clean_comment(item.comment);
        } else if (contains(item.decl, "using ")) {
            if (contains(item.decl, " = ")) {
                item.type = "Typedef";
                item.name = item.decl;
                item.name = partition(item.name, " = ")[0];
                item.name = replace(item.name, "using ", "");
                item.name = strip(item.name);
                item.comment = clean_comment(item.comment);
            } else {
                item.type = "Function Alias";
                item.name = partition(item.decl, "::")[2];
                item.name = replace(item.name, ";", "");
                item.name = strip(item.name);
                item.name += "()";
                item.comment = clean_comment(item.comment);
            }
        } else if (contains(item.decl, "enum ")) {
            item.type = "Enum";
            item.name = item.decl;
            item.name = replace(item.name, "enum ", "");
            item.name = replace(item.name, "struct ", "");
            item.name = replace(item.name, "{", "");
            item.name = strip(item.name);
            item.comment = clean_comment(item.comment);
            if (!item.children.empty()) {
                item.comment += "\n\n- Values:\n";
                for (auto& child : item.children) {
                    child.decl = replace(child.decl, ";", "");
                    child.decl = replace(child.decl, "}", "");
                    child.name = split(partition(child.decl, "=")[0]).back();
                    child.name = replace(child.name, ",", "");
                    item.comment += "    - " + child.name + ": " +
                                    replace(child.comment, "///", "");
                    item.decl += child.decl;
                }
                item.decl += "}\n";
            }
        } else if (contains(item.decl, "struct ")) {
            item.type = "Struct";
            item.name = item.decl;
            if (contains(item.name, "template ")) {
                item.name = partition(item.name, "\n")[2];
            }
            if (contains(item.name, " : ")) {
                item.name = partition(item.name, " : ")[0];
            }
            item.name = replace(item.name, "struct ", "");
            item.name = replace(item.name, "{", "");
            item.name = replace(item.name, ";", "");
            item.name = strip(item.name);
            item.comment = clean_comment(item.comment);
            if (!item.children.empty()) {
                item.comment += "\n\n- Members:\n";
                for (auto& child : item.children) {
                    auto isvar =
                        !contains(child.decl, " operator") &&
                        (!contains(child.decl, "(") ||
                            (child.decl.find("=") < child.decl.find("(")));
                    if (isvar) {
                        child.name =
                            replace(split(partition(child.decl, "=")[0]).back(),
                                ";", "");
                        item.comment += "    - " + child.name + ": " +
                                        replace(child.comment, "///", "");
                        item.decl += child.decl;
                    } else {
                        if (contains(child.decl, "{")) {
                            child.decl = partition(child.decl, "{")[0] + ";\n";
                        }
                        if (contains(child.decl, " : ")) {
                            child.decl =
                                partition(child.decl, " : ")[0] + ";\n";
                        }
                        child.decl = replace(child.decl, "\n", " ") + "\n";
                        while (contains(child.decl, "  "))
                            child.decl = replace(child.decl, "  ", " ");
                        child.decl = "   " + child.decl;
                        child.decl = replace(child.decl, " ;", ";");
                        child.name =
                            split(partition(child.decl, "(")[0]).back() + "()";
                        if (contains(child.decl, "operator "))
                            child.name = "operator " + child.name;
                        item.comment += "    - " + child.name + ": " +
                                        replace(child.comment, "///", "");
                        item.decl += child.decl;
                    }
                }
                item.decl += "}\n";
            }
        } else {
            auto isvar = !contains(item.decl, " operator") &&
                         (!contains(item.decl, "(") ||
                             (item.decl.find("=") < item.decl.find("(")));
            if (isvar) {
                if (contains(item.decl, "const ")) {
                    item.type = "Constant";
                } else {
                    item.type = "Variable";
                }
                item.name = split(partition(item.decl, "=")[0]).back();
                item.comment = clean_comment(item.comment);
            } else {
                item.type = "Function";
                if (contains(item.decl, "{")) {
                    item.decl = partition(item.decl, "{")[0] + ";\n";
                }
                if (contains(item.decl, " : ")) {
                    item.decl = partition(item.decl, " : ")[0] + ";\n";
                }
                // item.decl = replace_str(item.decl, "\n", " ") + "\n";
                // while(contains(item.decl, "  ")) item.decl =
                // replace_str(item.decl, "  ", " ");
                item.decl = replace(item.decl, "( ", "(");
                item.decl = replace(item.decl, " ;", ";");
                item.name = split(partition(item.decl, "(")[0]).back() + "()";
                if (contains(item.decl, "operator "))
                    item.name = "operator " + item.name;
                item.comment = clean_comment(item.comment);
            }
        }
    }

    // render comment blocks
    auto md = std::string();
    for (auto&& item : items) {
        if (item.name != "") {
            if (item.type == "Namespace") {
                md += "## ";
            } else {
                md += "### ";
            }
            md += item.type + " ";
            md += replace(replace(item.name, "<", " <"), ">", " \\>");
            // md += item.name;
            md += "\n\n";
        }
        if (item.decl != "") {
            md += "~~~ .cpp\n";
            md += replace(replace(item.decl, "{{", "{ {"), "}}", "} }");
            md += "~~~\n\n";
        }
        md += item.comment + "\n\n";
    }
    return md;
}

int main(int argc, char** argv) {
    for (auto filename : filenames) {
        auto cpp = load_txtfile(filename);
        auto md = make_doc(cpp);
        auto filename_out = filename;
        filename_out = replace(filename_out, ".h", ".md");
        filename_out = replace(filename_out, "yocto/", "docs/");
        save_txtfile(filename_out, md);
    }
}
