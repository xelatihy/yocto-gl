#! /usr/bin/env python3 -b

# extract markdown documentation from cpp
# complete hack

class Item:
    def __init__(self):
        self.name = ''
        self.decl = ''
        self.comment = ''
        self.type = ''
        self.children = []

def make_doc(cpp, first_only=False):
    # comment blocks
    items = []
    cur_item = None
    first = True
    indented = False
    enum_last = False
    for line in cpp.splitlines(True):
        if cur_item:
            if '///' in line:
                cur_item.comment += line
            elif first:
                cur_item = None
                first = False
            else:
                cur_item.decl += line
                if ';' in line or '{' in line:
                    cur_item = None
                if enum_last and ',' in line:
                    cur_item = None
        else:
            if not '///' in line: continue
            if line.startswith("    "):
                items[-1].children += [Item()]
                cur_item = items[-1].children[-1]
                cur_item.comment += line
                indented = True
                enum_last = 'enum ' in items[-1].decl
            else:
                items += [Item()]
                cur_item = items[-1]
                cur_item.comment += line
                indented = False
                enum_last = False

    def clean_comment(comment):
        if comment.startswith("    "):
            comment = comment.replace("    ///", "")
            comment = comment.replace("    /// ", "")
            comment = comment.strip()
        else:
            comment = comment.replace("/// ", "")
            comment = comment.replace("///", "")
            comment = comment.strip()
        return comment

    # main namespace
    main_namespace = ''

    # hack find
    def hack_find(str, s):
        return str.find(s) if str.find(s) >= 0 else 10000

    # fix type
    for item in items:
        if item.decl == "":
            item.comment = clean_comment(item.comment)
        elif "namespace " in item.decl:
            item.type = "Namespace"
            item.name = item.decl
            item.name = item.name.replace("namespace ", "").replace("{", "").strip()
            if (main_namespace == ""):
                main_namespace = item.name
            else:
                item.name = main_namespace + "::" + item.name
            item.decl = ""
            item.comment = clean_comment(item.comment)
        elif "using " in item.decl:
            if " = " in item.decl:
                item.type = "Typedef"
                item.name = item.decl
                item.name = item.name.partition(" = ")[0].replace("using ", "").strip()
                item.comment = clean_comment(item.comment)
            else:
                item.type = "Function Alias"
                item.name = item.decl.partition("::")[2].replace(";", "").strip() + "()"
                item.comment = clean_comment(item.comment)
        elif 'enum ' in item.decl:
            item.type = "Enum"
            item.name = item.decl.replace("enum ", "").replace("struct ", "").replace("{", "").strip()
            item.comment = clean_comment(item.comment)
            if item.children:
                item.comment += "\n\n- Values:\n"
                for child in item.children:
                    child.decl = child.decl.replace(";", "").replace("}", "")
                    child.name = child.decl.partition("=")[0].split()[-1].replace(",", "")
                    item.comment += "    - " + child.name + ": " + child.comment.replace("///", "")
                    item.decl += child.decl
                item.decl += "}\n"
        elif 'struct ' in item.decl:
            item.type = "Struct"
            item.name = item.decl
            if "template " in item.name:
                item.name = item.name.partition("\n")[2]
            if " : " in item.name:
                item.name = item.name.partition(" : ")[0]
            item.name = item.name.replace("struct ", "").replace("{", "").replace(";", "").strip()
            item.comment = clean_comment(item.comment)
            if item.children:
                item.comment += "\n\n- Members:\n"
                for child in item.children:
                    isvar = " operator" not in child.decl and ("(" not in child.decl or hack_find(child.decl, "=") < hack_find(child.decl, "("))
                    if isvar:
                        child.name = child.decl.partition("=")[0].split()[-1].replace(";", "")
                        item.comment += "    - " + child.name + ": " + child.comment.replace("///", "")
                        item.decl += child.decl
                    else:
                        if "{" in child.decl:
                            child.decl = child.decl.partition("{")[0] + ";\n"
                        if " : " in child.decl:
                            child.decl = child.decl.partition(" : ")[0] + ";\n"
                        child.decl = child.decl.replace("\n", " ") + "\n"
                        while '  ' in child.decl:
                            child.decl = child.decl.replace("  ", " ")
                        child.decl = "   " + child.decl.replace(" ;", ";")
                        child.name = child.decl.partition("(")[0].split()[-1] + "()"
                        if "operator " in child.decl:
                            child.name = "operator " + child.name
                        item.comment += "    - " + child.name + ": " + child.comment.replace("///", "")
                        item.decl += child.decl
                item.decl += "}\n"
        else:
            isvar = " operator" not in item.decl and ("(" not in item.decl or hack_find(item.decl, "=") < hack_find(item.decl, "("))
            if isvar:
                if  "const " in item.decl:
                    item.type = "Constant"
                else:
                    item.type = "Variable"
                item.name = item.decl.partition("=")[0].split()[-1]
                item.comment = clean_comment(item.comment)
            else:
                item.type = "Function"
                if "{" in item.decl:
                    item.decl = item.decl.partition("{")[0] + ";\n"
                if " : " in item.decl:
                    item.decl = item.decl.partition(" : ")[0] + ";\n"
                # item.decl = replace_str(item.decl, "\n", " ") + "\n"
                # while(contains(item.decl, "  ")) item.decl =
                # replace_str(item.decl, "  ", " ")
                item.decl = item.decl.replace("( ", "(").replace(" ;", ";")
                item.name = item.decl.partition("(")[0].split()[-1] + "()"
                if "operator " in item.decl:
                    item.name = "operator " + item.name
                item.comment = clean_comment(item.comment)

    # render comment blocks
    md = ''
    first = True
    for item in items:
        if item.name != "":
            md += "#### "
            md += item.type + " "
            md += item.name.replace("<", " <").replace(">", " \\>")
            # md += item.name
            md += "\n\n"
        if item.decl != "":
            md += "~~~ .cpp\n"
            md += item.decl.replace("{{", "{ {").replace("}}", "} }")
            md += "~~~\n\n"
        md += item.comment + "\n\n"
        if first:
            if first_only: return md
            md += "## API Documentation\n\n"
            first = False
    return md

template = '''
    <!DOCTYPE html>
    <html lang="en">
    <head>
      <title>Yocto/GL</title>
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=yes">
      <link rel="stylesheet" href="style.css">
    </head>
    <body>
    <header>
        <nav>
            <img src="images/logo_white.png">
            <a href="index.html">about</a>
            <a href="yocto_gl.html">api</a>
            <a href="https://github.com/xelatihy/yocto-gl">github</a>
        </nav>
    </header>
    <article>
    $body$
    <article>
    <footer></footer>
    </body>
    </html>
'''

def make_html(md):
    import markdown, glob
    html = markdown.markdown(md, ['markdown.extensions.extra',
                     'markdown.extensions.codehilite'],
                     output_format='html5')
    html = html.replace('<pre>', '<pre><code>')
    html = html.replace('</pre>', '</code></pre>')
    for link in glob.glob('docs/*.md'):
        link = link.replace('docs/','')
        hlink = link.replace('.md', '.html')
        html = html.replace(link, hlink)
    while '<p><img' in html:
        before, _, remainder = html.partition('<p><img')
        middle, _, after = remainder.partition('</p>')
        html = before + '<figure><img' + middle + '</figure>' + after
    html = template.replace('$body$', html)
    return html

for filename in ["yocto/yocto_gl.h", "yocto/yocto_gltf.h"]:
    with open(filename) as f: cpp = f.read()
    md = make_doc(cpp)
    html = make_html(md)
    filename_md = filename.replace(".h", ".md").replace("yocto/", "docs/")
    with open(filename_md, 'wt') as f: f.write(md)
    filename_html = filename_md.replace(".md", ".html")
    with open(filename_html, 'wt') as f: f.write(html)

for filename in ["yocto/yocto_gl.h"]:
    with open(filename) as f: cpp = f.read()
    md = make_doc(cpp, True)
    html = make_html(md)
    filename_md = filename.replace(".h", ".md").replace("yocto/", "docs/")
    with open('readme.md', 'wt') as f: f.write(md)
    with open('docs/index.html', 'wt') as f: f.write(html)
