# Yocto/Gui

A set of utilities for windows and widgets using GLFW and ImGui.

This library depends in yocto_math.h
The library depends on GLEW for OpenGL functions on Windows and Linux,
GLFW for window management and dear ImGui for widgets.


## History

- v 0.0: initial release split from yocto_glu

## Namespace ygui

Window Management and Immediate Mode Widgets

### Struct window

~~~ .cpp
struct window;
~~~

Forward declaration of window

### Function void()

~~~ .cpp
typedef void (*text_callback)(window*, unsigned int);
~~~

Text callback

### Function void()

~~~ .cpp
typedef void (*mouse_callback)(window*, int button, bool press, int mods);
~~~

Mouse callback

### Function void()

~~~ .cpp
typedef void (*refresh_callback)(window*);
~~~

Window refresh callback

### Function init_window()

~~~ .cpp
window* init_window(int width, int height, const std::string& title,
    void* user_pointer = nullptr);
~~~

initialize glfw

### Function clear_window()

~~~ .cpp
void clear_window(window* window);
~~~

Clear glfw

### Function set_callbacks()

~~~ .cpp
void set_callbacks(window* win, text_callback text_cb, mouse_callback mouse_cb,
    refresh_callback refresh_cb);
~~~

initialize glfw

### Function set_window_title()

~~~ .cpp
void set_window_title(window* win, const std::string& title);
~~~

Set window title

### Function wait_events()

~~~ .cpp
void wait_events(window* win);
~~~

Wait events

### Function poll_events()

~~~ .cpp
void poll_events(window* win);
~~~

Poll events

### Function swap_buffers()

~~~ .cpp
void swap_buffers(window* win);
~~~

Swap buffers

### Function should_close()

~~~ .cpp
bool should_close(window* win);
~~~

Should close

### Function get_user_pointer()

~~~ .cpp
void* get_user_pointer(window* window);
~~~

User pointer

### Function get_mouse_button()

~~~ .cpp
int get_mouse_button(window* window);
~~~

Mouse button

### Function get_mouse_posi()

~~~ .cpp
ym::vec2i get_mouse_posi(window* window);
~~~

Mouse position

### Function get_mouse_posf()

~~~ .cpp
ym::vec2f get_mouse_posf(window* window);
~~~

Mouse position

### Function get_window_size()

~~~ .cpp
ym::vec2i get_window_size(window* window);
~~~

Window size

### Function get_framebuffer_size()

~~~ .cpp
ym::vec2i get_framebuffer_size(window* window);
~~~

Framebuffer size

### Function get_screenshot()

~~~ .cpp
std::vector<ym::vec4b> get_screenshot(
    window* win, ym::vec2i& wh, bool flipy = true, bool back = false);
~~~

Read pixels

### Function init_widgets()

~~~ .cpp
void init_widgets(window* win);
~~~

Init ui

### Function clear_widgets()

~~~ .cpp
void clear_widgets(window* win);
~~~

Clear ui

### Function begin_widgets()

~~~ .cpp
bool begin_widgets(window* win, const std::string& title);
~~~

Begin draw widget

### Function end_widgets()

~~~ .cpp
void end_widgets(window* win);
~~~

End draw widget

### Function get_widget_active()

~~~ .cpp
bool get_widget_active(window* win);
~~~

Whether widget are active

### Function separator_widget()

~~~ .cpp
void separator_widget(window* win);
~~~

Horizontal separator

### Function indent_begin_widgets()

~~~ .cpp
void indent_begin_widgets(window* win);
~~~

Indent widget

### Function indent_end_widgets()

~~~ .cpp
void indent_end_widgets(window* win);
~~~

Indent widget

### Function label_widget()

~~~ .cpp
void label_widget(window* win, const std::string& lbl, const std::string& mesg);
~~~

Label widget

### Function label_widget()

~~~ .cpp
void label_widget(
    window* win, const std::string& lbl, int val, const char* fmt = "%d");
~~~

Label widget

### Function label_widget()

~~~ .cpp
void label_widget(window* win, const std::string& lbl, ym::vec2i val,
    const char* fmt = "[ %d, %d ]");
~~~

Label widget

### Function label_widget()

~~~ .cpp
void label_widget(window* win, const std::string& lbl, ym::vec3i val,
    const char* fmt = "[ %d, %d, %d ]");
~~~

Label widget

### Function label_widget()

~~~ .cpp
void label_widget(window* win, const std::string& lbl, ym::vec4i val,
    const char* fmt = "[ %d, %d, %d, %d ]");
~~~

Label widget

### Function label_widget()

~~~ .cpp
void label_widget(
    window* win, const std::string& lbl, float val, const char* fmt = "%f");
~~~

Label and float widget

### Function label_widget()

~~~ .cpp
void label_widget(window* win, const std::string& lbl, ym::vec2f val,
    const char* fmt = "[ %f, %f ]");
~~~

Label and float widget

### Function label_widget()

~~~ .cpp
void label_widget(window* win, const std::string& lbl, ym::vec3f val,
    const char* fmt = "[ %f, %f, %f ]");
~~~

Label and float widget

### Function label_widget()

~~~ .cpp
void label_widget(window* win, const std::string& lbl, ym::vec4f val,
    const char* fmt = "[ %f, %f, %f ]");
~~~

Label and float widget

### Function slider_widget()

~~~ .cpp
bool slider_widget(window* win, const std::string& lbl, int* val, int min,
    int max, int incr = 1);
~~~

Int widget

### Function slider_widget()

~~~ .cpp
bool slider_widget(window* win, const std::string& lbl, ym::vec2i* val, int min,
    int max, int incr = 1);
~~~

Int widget

### Function slider_widget()

~~~ .cpp
bool slider_widget(window* win, const std::string& lbl, ym::vec3i* val, int min,
    int max, int incr = 1);
~~~

Int widget

### Function slider_widget()

~~~ .cpp
bool slider_widget(window* win, const std::string& lbl, ym::vec4i* val, int min,
    int max, int incr = 1);
~~~

Int widget

### Function slider_widget()

~~~ .cpp
bool slider_widget(window* win, const std::string& lbl, float* val, float min,
    float max, float incr = 1.0f);
~~~

Float widget

### Function slider_widget()

~~~ .cpp
bool slider_widget(window* win, const std::string& lbl, ym::vec2f* val,
    float min, float max, float incr = 1.0f);
~~~

Float widget

### Function slider_widget()

~~~ .cpp
bool slider_widget(window* win, const std::string& lbl, ym::vec3f* val,
    float min, float max, float incr = 1.0f);
~~~

Float widget

### Function slider_widget()

~~~ .cpp
bool slider_widget(window* win, const std::string& lbl, ym::vec4f* val,
    float min, float max, float incr = 1.0f);
~~~

Float widget

### Function checkbox_widget()

~~~ .cpp
bool checkbox_widget(window* win, const std::string& lbl, bool* val);
~~~

Bool widget

### Function combo_widget()

~~~ .cpp
bool combo_widget(window* win, const std::string& lbl, int* val,
    const std::vector<std::pair<std::string, int>>& labels);
~~~

Enum widget

### Function combo_widget()

~~~ .cpp
bool combo_widget(window* win, const std::string& lbl, void** val,
    const std::vector<std::pair<std::string, void*>>& labels);
~~~

Enum widget

### Function combo_widget()

~~~ .cpp
template <typename T>
inline bool combo_widget(window* win, const std::string& lbl, T* val,
    const std::vector<std::pair<std::string, T>>& labels);
~~~

Enum widget (assume we can convert T to int)

### Function combo_widget()

~~~ .cpp
template <typename T>
inline bool combo_widget(window* win, const std::string& lbl, T** val,
    const std::vector<std::pair<std::string, T*>>& labels);
~~~

Enum widget

### Function list_widget()

~~~ .cpp
bool list_widget(window* win, const std::string& lbl, int* val,
    const std::vector<std::pair<std::string, int>>& labels);
~~~

List widget

### Function list_widget()

~~~ .cpp
bool list_widget(window* win, const std::string& lbl, void** val,
    const std::vector<std::pair<std::string, void*>>& labels);
~~~

List widget

### Function list_widget()

~~~ .cpp
template <typename T>
inline bool list_widget(window* win, const std::string& lbl, T** val,
    const std::vector<std::pair<std::string, T*>>& labels);
~~~

List widget

### Function button_widget()

~~~ .cpp
bool button_widget(window* win, const std::string& lbl);
~~~

Button widget

### Function collapsing_header_widget()

~~~ .cpp
bool collapsing_header_widget(window* win, const std::string& lbl);
~~~

Collapsible header

### Function tree_begin_widget()

~~~ .cpp
bool tree_begin_widget(window* win, const std::string& lbl);
~~~

Start tree node

### Function tree_end_widget()

~~~ .cpp
void tree_end_widget(window* win);
~~~

End tree widget

### Function tree_begin_widget()

~~~ .cpp
bool tree_begin_widget(
    window* win, const std::string& lbl, void** selection, void* content);
~~~

Start selectable tree node

### Function tree_end_widget()

~~~ .cpp
void tree_end_widget(window* win, void* content);
~~~

End selectable tree node

### Function tree_leaf_widget()

~~~ .cpp
void tree_leaf_widget(
    window* win, const std::string& lbl, void** selection, void* content);
~~~

Selectable tree leaf node

### Function image_widget()

~~~ .cpp
void image_widget(window* win, int tid, ym::vec2i size);
~~~

Image widget

### Function scroll_region_begin_widget()

~~~ .cpp
void scroll_region_begin_widget(
    window* win, const std::string& lbl, int height, bool border);
~~~

Scroll region

### Function scroll_region_end_widget()

~~~ .cpp
void scroll_region_end_widget(window* win);
~~~

Scroll region

### Function scroll_region_here_widget()

~~~ .cpp
void scroll_region_here_widget(window* win);
~~~

Scroll region

