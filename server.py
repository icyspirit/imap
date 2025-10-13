
from typing import Optional
from fastapi import FastAPI, Response, Query
from pydantic import BaseModel
import mapping


app = FastAPI()


class InputData(BaseModel):
    shot: int
    tbegin: Optional[float] = None
    tend: Optional[float] = None
    dd_version: str = "3.39.0"


@app.post("/{ids_name}")
def mapping_endpoint(
    ids_name: str,
    inp: InputData,
    postprocess: bool = Query(False, description="Postprocess after mapping"),
    tree: str = Query("", description="Mapping tree"),
):
    ids_url = ids_name

    query_parts = []
    if postprocess:
        query_parts.append("postprocess=true")
    if tree:
        query_parts.append(f"tree={tree}")
    if query_parts:
        ids_url += "?" + "&".join(query_parts)

    ids = mapping.generate_mapping(
        ids_url,
        inp.shot,
        inp.tbegin,
        inp.tend,
        inp.dd_version,
    )

    return Response(
        content=ids.serialize(),
        media_type="application/octet-stream",
    )


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app)
