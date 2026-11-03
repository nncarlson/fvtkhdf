"""Render method signatures above an indented reStructuredText description."""

from docutils import nodes
from sphinx.util.docutils import SphinxDirective


class MethodDoc(SphinxDirective):
    has_content = True

    def run(self):
        self.assert_has_content()
        try:
            separator = self.content.index("")
        except ValueError:
            raise self.error("Separate the signatures and description with a blank line.")
        if not separator or not any(self.content[separator + 1:]):
            raise self.error("Both signatures and a description are required.")

        text = "\n".join(self.content[:separator])
        signatures = nodes.literal_block(text, text, language="text")
        self.set_source_info(signatures)
        description = nodes.container(classes=["method-description"])
        self.state.nested_parse(
            self.content[separator + 1:],
            self.content_offset + separator + 1,
            description,
        )
        return [nodes.container("", signatures, description, classes=["method-doc"])]


def setup(app):
    app.add_directive("method-doc", MethodDoc)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
