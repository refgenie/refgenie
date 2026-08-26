from typing import Any

import jinja2

from refgenie.logger import logger


def jinja_render_template_strictly(template: str, values: Any) -> str:
    """
    Render a command template with StrictUndefined: every attribute the
    template references must exist on ``values``, otherwise ValueError is
    raised.

    Args:
        template: Jinja2 template string.
        values: The object exposed to the template as ``values``; the template
            reaches into it by attribute or item access.

    Returns:
        The rendered command string.
    """
    env = jinja2.Environment(
        undefined=jinja2.StrictUndefined,
    )
    templ_obj = env.from_string(template)
    logger.debug(f"Rendering template with values: {values}")
    try:
        rendered = templ_obj.render(dict(values=values))
    except jinja2.exceptions.UndefinedError as e:
        logger.debug(f"Template: '{template}'")
        raise ValueError(f"Error populating command template: {str(e)}") from e
    logger.debug(f"Rendered template: {rendered}")
    return rendered
