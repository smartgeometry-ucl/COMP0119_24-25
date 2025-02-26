"use strict";

class TutteEmbedding {

	constructor(geometry, boundaryscheme, weightset) {
		this.geometry = geometry;
		this.vertexNumber = geometry.mesh.vertices.length;

		this.computeLongestBoundary();
		this.computeBoundaryEdgeLength()
		this.computeUbarVbar(boundaryscheme);
		this.computeMatrix(weightset);
		this.computeUV();
	}

	computeLongestBoundary() {
		this.boundary = undefined;
		this.boundaryVertexNumber = 0;

		let lengthTmp = 0;
		for (let b of this.geometry.mesh.boundaries) {
			let length = 0
			for (let bv of b.adjacentVertices()) {
				length++;
			}
			if (length > lengthTmp) {
				this.boundary = b;
				this.boundaryVertexNumber = length;
			}
			lengthTmp = b.length;
		}
	}

	computeBoundaryEdgeLength() {
		this.boundaryEdgeLength = {};
		this.boundaryTotalLength = 0;

		for (let bhe of this.boundary.adjacentHalfedges()) {
			this.boundaryEdgeLength[bhe.vertex] = this.geometry.length(bhe.edge);
			this.boundaryTotalLength = this.boundaryTotalLength + this.boundaryEdgeLength[bhe.vertex];
		}
	}

	computeUbarVbar(boundaryscheme) {
		this.ubar = DenseMatrix.zeros(this.vertexNumber);
		this.vbar = DenseMatrix.zeros(this.vertexNumber);
		
		let radius = 1;
		//let radius = this.boundaryTotalLength / (2 * Math.PI);
		let l = 0;

		let u0 = 0;
		let v0 = 0;

		switch (boundaryscheme) {
			case "Even": default:
				let count = 0;
				for (let bv of this.boundary.adjacentVertices()) {
					let i = bv.index;
					u0 = radius * Math.cos(l/radius);
					v0 = -radius * Math.sin(l/radius);
					this.ubar.set(u0, i);
					this.vbar.set(v0, i);
					count++;
					l = count/this.boundaryVertexNumber * 2 * Math.PI * radius;
				}
				break;
			case "By Scale":
				let s = 0;
				for (let bv of this.boundary.adjacentVertices()) {
					let i = bv.index;
					
					
					const theta = l/radius;
					var sf = 1;
					const t = 1.0;
					
					
					// Please edit here to see if you can make the parametrisation map to a square instead of a circle.
					
					u0 = radius * Math.cos(theta) * sf;
					v0 = -radius * Math.sin(theta) * sf;
					
					
					
					this.ubar.set(u0, i);
					this.vbar.set(v0, i);
					s = s + this.boundaryEdgeLength[bv];
					l = s * 2 * Math.PI * radius / this.boundaryTotalLength;
				}
				break;
		}
	}

	computeMatrix(weightset) {
		let T = new Triplet(this.vertexNumber, this.vertexNumber);
		switch (weightset) {
			case "Uniform Laplacian": default:
				for (let v of this.geometry.mesh.vertices) {
					let i = v.index;
					if (this.ubar.get(i) != 0 || this.vbar.get(i) != 0) {
						T.addEntry(1, i, i);
					} else {
						let n = 0;
						for (let nv of v.adjacentVertices()) {
							let j = nv.index;
							let weight = 1;
							T.addEntry(weight, i, j);
							n = n + weight;
						}
						T.addEntry(-n, i, i);
					}
				}
				break;
			case "Laplace-Beltrami":
				for (let v of this.geometry.mesh.vertices) {
					let i = v.index;
					if (this.ubar.get(i) != 0 || this.vbar.get(i) != 0) {
						T.addEntry(1, i, i);
					} else {
						let n = 0;
						let area = this.geometry.circumcentricDualArea(v);
						for (let nhe of v.adjacentHalfedges()) {
							let j = nhe.next.vertex.index;
							let weight = (this.geometry.cotan(nhe) + this.geometry.cotan(nhe.twin)) / (2 * area);
							T.addEntry(weight, i, j);
							n = n + weight;
						}
						T.addEntry(-n, i, i);
					}
				}
				break;
			case "Mean Value":
				for (let v of this.geometry.mesh.vertices) {
					let i = v.index;
					if (this.ubar.get(i) != 0 || this.vbar.get(i) != 0) {
						T.addEntry(1, i, i);
					} else {
						let n = 0;
						for (let nhe of v.adjacentHalfedges()) {
							let j = nhe.next.vertex.index;
							let distance = this.geometry.length(nhe.edge);
							let delta = nhe.next.corner;
							let gamma = nhe.twin.prev.corner;
							let tdelta = delta ? Math.tan(this.geometry.angle(delta) / 2) : 0;
							let tgamma = gamma ? Math.tan(this.geometry.angle(gamma) / 2) : 0;
							let weight = (tdelta + tgamma) / distance;
							T.addEntry(weight, i, j);
							n = n + weight;
						}
						T.addEntry(-n, i, i);
					}
				}
				break;
		}

		this.A = SparseMatrix.fromTriplet(T);
	}

	computeUV() {
		let lu = this.A.lu();
		this.u = lu.solveSquare(this.ubar);
		this.v = lu.solveSquare(this.vbar);
	}

	apply() {
		let parameterization = {};
		for (let v of this.geometry.mesh.vertices) {
			parameterization[v] = new Vector(this.u.get(v.index), this.v.get(v.index));
		}

		//normalize(parameterization, this.geometry.mesh.vertices);
		
		return parameterization;
	}

}
